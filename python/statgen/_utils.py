import warnings

import numpy as np

CANONICAL_CHR_ORDER = [str(i) for i in range(1, 23)] + ["X"]
CHR_RANK = {c: i for i, c in enumerate(CANONICAL_CHR_ORDER)}
_ALLELE_HASH_P = 2_147_483_647
_ALLELE_HASH_BASE1 = 257
_ALLELE_HASH_BASE2 = 263
_ALLELE_HASH_MAX_CHARS = 150
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


def allele_hash64(alleles) -> np.ndarray:
    arr = np.asarray(alleles, dtype=object).reshape(-1)
    n = arr.size
    if n == 0:
        return np.array([], dtype=np.uint64)

    strings = [str(allele) for allele in arr]
    lengths = np.fromiter((len(s) for s in strings), dtype=np.int64, count=n)
    if np.any(lengths > _ALLELE_HASH_MAX_CHARS):
        warnings.warn(
            "Allele length exceeds 150 characters; hashing uses first 150 characters",
            RuntimeWarning,
            stacklevel=2,
        )
        strings = [s[:_ALLELE_HASH_MAX_CHARS] for s in strings]
        lengths = np.minimum(lengths, _ALLELE_HASH_MAX_CHARS)

    out = np.empty(n, dtype=np.uint64)
    single_base = lengths == 1
    if np.any(single_base):
        single_idx = np.flatnonzero(single_base)
        single_strings = np.array([strings[i] for i in single_idx], dtype=object)
        assigned = np.zeros(single_idx.size, dtype=bool)
        for base, hash_value in zip(("A", "C", "G", "T"), _SINGLE_BASE_HASH64):
            matches = single_strings == base
            if np.any(matches):
                out[single_idx[matches]] = hash_value
                assigned[matches] = True

        if np.any(~assigned):
            other_idx = single_idx[~assigned]
            encoded = [strings[i].encode("utf-8") for i in other_idx]
            h1, h2 = _allele_hash_lanes_from_utf8(encoded)
            out[other_idx] = (
                h1.astype(np.uint64) << np.uint64(32)
            ) + h2.astype(np.uint64)

    if np.all(single_base):
        return out

    multi_idx = np.flatnonzero(~single_base)
    encoded = [strings[i].encode("utf-8") for i in multi_idx]
    h1, h2 = _allele_hash_lanes_from_utf8(encoded)
    out[multi_idx] = (h1.astype(np.uint64) << np.uint64(32)) + h2.astype(np.uint64)
    return out


def _allele_hash_lanes_from_utf8(encoded) -> tuple[np.ndarray, np.ndarray]:
    n = len(encoded)
    byte_lengths = np.fromiter((len(b) for b in encoded), dtype=np.int64, count=n)
    max_len = int(byte_lengths.max(initial=0))
    mat = np.zeros((n, max_len), dtype=np.uint8)
    if max_len:
        flat = np.frombuffer(b"".join(encoded), dtype=np.uint8)
        row_idx = np.repeat(np.arange(n, dtype=np.int64), byte_lengths)
        starts = np.cumsum(byte_lengths, dtype=np.int64) - byte_lengths
        col_idx = np.arange(flat.size, dtype=np.int64) - np.repeat(starts, byte_lengths)
        mat[row_idx, col_idx] = flat

    h1 = np.ones(n, dtype=np.int64)
    h2 = np.ones(n, dtype=np.int64)
    for j in range(max_len):
        active = j < byte_lengths
        x = mat[active, j].astype(np.int64) + 1
        h1[active] = (h1[active] * _ALLELE_HASH_BASE1 + x) % _ALLELE_HASH_P
        h2[active] = (h2[active] * _ALLELE_HASH_BASE2 + x) % _ALLELE_HASH_P

    return h1, h2


_A_HASH64, _C_HASH64, _G_HASH64, _T_HASH64 = _SINGLE_BASE_HASH64


def _is_single_base(h) -> np.ndarray:
    return np.isin(h, _SINGLE_BASE_HASH64)


def _is_single_nucleotide_variant(a1_hash64, a2_hash64) -> np.ndarray:
    return np.isin(a1_hash64, _SINGLE_BASE_HASH64) & np.isin(
        a2_hash64, _SINGLE_BASE_HASH64
    )


def _is_strand_ambiguous(a1_hash64, a2_hash64) -> np.ndarray:
    a1 = np.asarray(a1_hash64, dtype=np.uint64).reshape(-1)
    a2 = np.asarray(a2_hash64, dtype=np.uint64).reshape(-1)
    return (
        ((a1 == _A_HASH64) & (a2 == _T_HASH64))
        | ((a1 == _T_HASH64) & (a2 == _A_HASH64))
        | ((a1 == _C_HASH64) & (a2 == _G_HASH64))
        | ((a1 == _G_HASH64) & (a2 == _C_HASH64))
    )
