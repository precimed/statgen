import warnings

import numpy as np

CANONICAL_CHR_ORDER = [str(i) for i in range(1, 23)] + ["X"]
CHR_RANK = {c: i for i, c in enumerate(CANONICAL_CHR_ORDER)}
_ALLELE_HASH_P = 2_147_483_647
_ALLELE_HASH_BASE1 = 257
_ALLELE_HASH_BASE2 = 263
_ALLELE_HASH_MAX_CHARS = 1_000_000
# Includes flat bytes, int64 powers/products/indices and reduction buffers.
_ALLELE_HASH_BATCH_BYTES = 262_144
_SINGLE_BASE_BYTES = np.frombuffer(b"ACGT", dtype=np.uint8).astype(np.int64) + 1
_SINGLE_BASE_HASH64 = (
    (
        ((_ALLELE_HASH_BASE1 + _SINGLE_BASE_BYTES) % _ALLELE_HASH_P).astype(np.uint64)
        << np.uint64(32)
    )
    + ((_ALLELE_HASH_BASE2 + _SINGLE_BASE_BYTES) % _ALLELE_HASH_P).astype(np.uint64)
)


def validate_requested_shards(shards, available_labels, where: str) -> list[str]:
    labels = list(available_labels)
    if shards is None:
        return labels
    if isinstance(shards, str):
        raise ValueError(f"{where}: shards must be a non-empty list of unique canonical contig labels")

    requested = list(shards)
    if not requested:
        raise ValueError(f"{where}: shards must be a non-empty list of unique canonical contig labels")

    seen = set()
    prev_rank = None
    for label in requested:
        if label not in CHR_RANK:
            raise ValueError(f"{where}: unsupported shard label {label!r}; expected canonical labels 1-22 or X")
        if label in seen:
            raise ValueError(f"{where}: duplicate shard label {label!r} in shards")
        rank = CHR_RANK[label]
        if prev_rank is not None and rank <= prev_rank:
            raise ValueError(f"{where}: shards must be in canonical subsequence order")
        seen.add(label)
        prev_rank = rank
        if label not in labels:
            raise ValueError(f"{where}: requested shard {label!r} is not present")
    return requested


def _warn_allele_cap(a1, a2, path, rows):
    """Source-boundary diagnostic; helpers deliberately remain silent."""
    lengths = (a1.str.len().to_numpy(), a2.str.len().to_numpy())
    affected = tuple(lens > _ALLELE_HASH_MAX_CHARS for lens in lengths)
    count = sum(np.count_nonzero(mask) for mask in affected)
    if count:
        first = min((int(np.flatnonzero(mask)[0]), col)
                    for col, mask in enumerate(affected) if mask.any())
        i, col = first
        warnings.warn(
            f"{path}: {count} allele(s) exceed hash cap {_ALLELE_HASH_MAX_CHARS}; "
            f"first affected row {rows[i]}, column a{col + 1}, length {lengths[col][i]}; "
            f"hashing uses first {_ALLELE_HASH_MAX_CHARS} bases",
            RuntimeWarning, stacklevel=2,
        )


def _allele_hash_powers(base, length):
    powers = np.empty(length + 1, dtype=np.int64)
    powers[0] = 1
    used = 1
    factor = base
    # PERF: double the initialized prefix with native integer operations.
    # Products of residues are below (2^31 - 1)^2 < 2^62.
    while used < powers.size:
        count = min(used, powers.size - used)
        powers[used:used + count] = powers[:count] * factor % _ALLELE_HASH_P
        factor = factor * factor % _ALLELE_HASH_P
        used += count
    return powers


def allele_hash64(alleles) -> np.ndarray:
    arr = np.asarray(alleles, dtype=object).reshape(-1)
    if not arr.size:
        return np.empty(0, dtype=np.uint64)
    # Preserve lazy pandas loading for LD/cache-only workflows.
    import pandas as pd

    strings = pd.Series(arr, copy=False)
    mapped = strings.map(dict(zip(("A", "C", "G", "T"), _SINGLE_BASE_HASH64)))
    idx = np.flatnonzero(mapped.isna().to_numpy())
    # Single-base hashes are below 2^53, so map's float64 missing-value
    # representation preserves every bit before conversion to uint64.
    out = mapped.to_numpy(dtype=np.uint64, na_value=0)
    if not idx.size:
        return out
    strings = strings.iloc[idx].reset_index(drop=True)
    lengths = strings.str.len().to_numpy(dtype=np.int64)
    lengths = np.minimum(lengths, _ALLELE_HASH_MAX_CHARS)
    longest = int(lengths.max(initial=0))
    powers = [_allele_hash_powers(b, longest)
              for b in (_ALLELE_HASH_BASE1, _ALLELE_HASH_BASE2)]
    ends = np.cumsum(lengths, dtype=np.int64)
    # Native linear partitioning. Grouping by end offsets bounds small batches
    # by twice the target; larger individual alleles occupy their own batch.
    target = max(1, _ALLELE_HASH_BATCH_BYTES // 2)
    large = lengths > target
    split = np.r_[True, np.diff((ends - 1) // target) != 0] | large | np.r_[False, large[:-1]]
    batch_starts = np.flatnonzero(split)
    batch_ends = np.r_[batch_starts[1:], idx.size]
    # PERF: batches bound all byte expansions; isolated long alleles use blocks.
    for start, stop in zip(batch_starts, batch_ends):
        lens = lengths[start:stop]
        if lens[0] > _ALLELE_HASH_BATCH_BYTES:
            allele = strings.iat[start]
            lanes = [int(power[lens[0]]) for power in powers]
            for offset in range(0, int(lens[0]), _ALLELE_HASH_BATCH_BYTES):
                piece = allele[offset:min(offset + _ALLELE_HASH_BATCH_BYTES, int(lens[0]))]
                x = np.frombuffer(piece.encode("ascii"), dtype=np.uint8).astype(np.int64) + 1
                positions = lens[0] - offset - 1 - np.arange(x.size)
                for lane, power in enumerate(powers):
                    lanes[lane] = (lanes[lane] + int(np.sum(x * power[positions], dtype=np.int64))) % _ALLELE_HASH_P
            out[idx[start]] = np.uint64((lanes[0] << 32) + lanes[1])
        else:
            # Slice before encoding; valid DNA has one byte per character.
            encoded = strings.iloc[start:stop].str.slice(0, _ALLELE_HASH_MAX_CHARS).str.encode("ascii")
            x = np.frombuffer(b"".join(encoded), dtype=np.uint8).astype(np.int64) + 1
            starts = np.cumsum(lens) - lens
            positions = np.repeat(starts + lens, lens) - np.arange(x.size) - 1
            lanes = [((np.add.reduceat(x * power[positions], starts) + power[lens])
                      % _ALLELE_HASH_P).astype(np.uint64) for power in powers]
            out[idx[start:stop]] = (lanes[0] << np.uint64(32)) + lanes[1]
    return out


_A_HASH64, _C_HASH64, _G_HASH64, _T_HASH64 = _SINGLE_BASE_HASH64


def _is_single_base(h) -> np.ndarray:
    return np.isin(h, _SINGLE_BASE_HASH64)


def _is_single_nucleotide_variant(a1_hash64, a2_hash64) -> np.ndarray:
    return _is_single_base(a1_hash64) & _is_single_base(a2_hash64)


def _is_strand_ambiguous(a1_hash64, a2_hash64) -> np.ndarray:
    a1 = np.asarray(a1_hash64, dtype=np.uint64).reshape(-1)
    a2 = np.asarray(a2_hash64, dtype=np.uint64).reshape(-1)
    return (
        ((a1 == _A_HASH64) & (a2 == _T_HASH64))
        | ((a1 == _T_HASH64) & (a2 == _A_HASH64))
        | ((a1 == _C_HASH64) & (a2 == _G_HASH64))
        | ((a1 == _G_HASH64) & (a2 == _C_HASH64))
    )


def _validate_allow_indel_swaps(value) -> bool:
    if not isinstance(value, (bool, np.bool_)):
        raise ValueError("allow_indel_swaps must be a boolean scalar")
    return bool(value)
