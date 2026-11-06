from dataclasses import dataclass

import numpy as np

from ._utils import _is_single_base, _is_single_nucleotide_variant, _is_strand_ambiguous

_KEY_DTYPE = np.dtype([("bp", "<i8"), ("a1_hash64", "<u8"), ("a2_hash64", "<u8")])


def numeric_variant_keys(bp_arr, a1_hash64, a2_hash64) -> np.ndarray:
    bp = np.asarray(bp_arr, dtype=np.int64).reshape(-1)
    a1h = np.asarray(a1_hash64, dtype=np.uint64).reshape(-1)
    a2h = np.asarray(a2_hash64, dtype=np.uint64).reshape(-1)
    if not (bp.size == a1h.size == a2h.size):
        raise ValueError("numeric variant key vector lengths mismatch")
    keys = np.empty(bp.size, dtype=_KEY_DTYPE)
    keys["bp"] = bp
    keys["a1_hash64"] = a1h
    keys["a2_hash64"] = a2h
    return keys


def check_unique_sorted_variant_keys(sorted_keys: np.ndarray, where: str, object_name: str = "source") -> None:
    if sorted_keys.size < 2:
        return
    dup = sorted_keys[1:] == sorted_keys[:-1]
    if np.any(dup):
        pos = int(np.flatnonzero(dup)[0] + 1)
        key = sorted_keys[pos]
        raise ValueError(
            f"Ambiguous duplicate {object_name}/reference matching key in {where}: "
            f"bp={int(key['bp'])}, a1_hash64={int(key['a1_hash64'])}, "
            f"a2_hash64={int(key['a2_hash64'])}"
        )


@dataclass(frozen=True)
class VariantMatch:
    source_rows: np.ndarray
    is_allele_swapped: np.ndarray
    num_strand_ambiguous: int
    num_disabled_indels: int


def _lookup_candidates(ref_keys, src_keys, source_rows, shard_label, source_name, candidate_eligible=None):
    """Return unique eligible source rows and count hits from ineligible rows."""
    order = np.argsort(src_keys, kind="mergesort", order=("bp", "a1_hash64", "a2_hash64"))
    sorted_keys = src_keys[order]
    left = np.searchsorted(sorted_keys, ref_keys, side="left")
    right = np.searchsorted(sorted_keys, ref_keys, side="right")
    present = right > left
    accepted = present.copy()
    if candidate_eligible is not None:
        # A key fixes both allele hashes, so all candidates for a slot have
        # the same allele class and eligibility, even when they are duplicates.
        accepted[present] = candidate_eligible[order[left[present]]]
    excluded_hits = int(np.sum((right - left)[present & ~accepted]))
    ambiguous = (right - left > 1) & accepted
    if np.any(ambiguous):
        key = ref_keys[np.flatnonzero(ambiguous)[0]]
        raise ValueError(
            f"Ambiguous duplicate {source_name}/reference matching key in shard {shard_label}: "
            f"bp={int(key['bp'])}, a1_hash64={int(key['a1_hash64'])}, "
            f"a2_hash64={int(key['a2_hash64'])}"
        )
    out = np.full(ref_keys.size, -1, dtype=np.int64)
    out[accepted] = source_rows[order[left[accepted]]]
    return out, excluded_hits


def match_shard_numeric(
    ref_keys: np.ndarray,
    src_keys: np.ndarray,
    shard_label: str,
    *,
    source_name: str = "source",
    allow_indel_swaps: bool = False,
) -> VariantMatch:
    # Reference validation belongs to load_reference; cached references are trusted.
    rows, _ = _lookup_candidates(
        ref_keys, src_keys, np.arange(src_keys.size), shard_label, source_name
    )
    swapped = np.zeros(ref_keys.size, dtype=bool)
    unused = np.ones(src_keys.size, dtype=bool)
    unused[rows[rows >= 0]] = False
    a1h = src_keys["a1_hash64"]
    a2h = src_keys["a2_hash64"]
    snv = _is_single_nucleotide_variant(a1h, a2h)
    indel = _is_single_base(a1h) ^ _is_single_base(a2h)
    eligible = snv | (indel & allow_indel_swaps)
    source_rows = np.flatnonzero(unused & (snv | indel))
    reversed_keys = numeric_variant_keys(
        src_keys["bp"][source_rows], a2h[source_rows], a1h[source_rows]
    )
    missing = np.flatnonzero(rows < 0)
    fallback, disabled = _lookup_candidates(
        ref_keys[missing], reversed_keys, source_rows, shard_label, source_name,
        candidate_eligible=eligible[source_rows],
    )
    accepted = missing[fallback >= 0]
    rows[accepted] = fallback[fallback >= 0]
    swapped[accepted] = True
    selected = rows[swapped]
    ambiguous = _is_strand_ambiguous(a1h[selected], a2h[selected])
    return VariantMatch(rows, swapped, int(np.count_nonzero(ambiguous)), disabled)
