"""Exact ID coding and numeric allele-group matching shared by source loaders."""

import warnings

import numpy as np
import pandas as pd

from ._utils import _validate_allow_indel_swaps, allele_hash64, _is_single_nucleotide_variant, _is_strand_ambiguous
from ._variant_match import numeric_variant_keys
from ._verbosity import info


def validate_matching_options(match_by, allow_indel_swaps):
    allow_indel_swaps = _validate_allow_indel_swaps(allow_indel_swaps)
    if not isinstance(match_by, str) or match_by not in ("genomic_position", "id"):
        raise ValueError('match_by must be "genomic_position" or "id"')
    if match_by == "id" and allow_indel_swaps:
        raise ValueError('allow_indel_swaps=true is incompatible with match_by="id"')
    return allow_indel_swaps


def _eligible(ids, a1, a2, context):
    snv = _is_single_nucleotide_variant(a1, a2)
    palindromic = _is_strand_ambiguous(a1, a2)
    for mask, label in ((palindromic, "palindromic SNV"), (~snv, "multi-base variant")):
        if mask.any():
            info(f"{context}: excluded {np.count_nonzero(mask)} {label}(s)")
    alleles_ok = snv & ~palindromic & (a1 != a2)
    names = pd.Series(ids, dtype=object)
    missing = (names.isna() | names.isin(["", "."])).to_numpy()
    valid = names.str.fullmatch(r"[\x21-\x7e]{1,32}", na=False).to_numpy(dtype=bool) & ~missing
    for mask, label in ((missing, "missing/empty/."), (~valid & ~missing, "overlength or invalid-character")):
        count = np.count_nonzero(mask & alleles_ok)
        if count:
            warnings.warn(f"{context}: excluded {count} unusable {label} ID(s)",
                          RuntimeWarning, stacklevel=2)
    return alleles_ok & valid


def match_id_panel(reference, sources, *, chromosome_scope=False):
    """Return panel-axis source indices and inversion flags.

    sources contains (context, IDs, a1 hashes, a2 hashes, chromosome labels)
    tuples. Returned indices address the concatenation of those source arrays.
    """
    ref_a1 = reference.a1_hash64
    ref_a2 = reference.a2_hash64
    ref_ids = reference.snp
    ref_ok = _eligible(ref_ids, ref_a1, ref_a2, "reference panel")
    if chromosome_scope:
        labels = pd.Index([s.label for s in reference.shards])
    source_ok = []
    for context, ids, a1, a2, chromosomes in sources:
        if chromosome_scope:
            scope = labels.get_indexer(chromosomes) >= 0
            eligible = np.zeros(len(ids), dtype=bool)
            eligible[scope] = _eligible(ids[scope], a1[scope], a2[scope], context)
        else:
            eligible = _eligible(ids, a1, a2, context)
        source_ok.append(eligible)
    src_ids = np.concatenate([s[1] for s in sources])
    src_a1 = np.concatenate([s[2] for s in sources])
    src_a2 = np.concatenate([s[3] for s in sources])
    src_ok = np.concatenate(source_ok)
    ri = np.flatnonzero(ref_ok)
    si = np.flatnonzero(src_ok)
    codes, _ = pd.factorize(np.concatenate((ref_ids[ri], src_ids[si])), sort=False)
    # A,C,G,T -> 0,1,2,3. Complement is 3-code.
    hashes = allele_hash64(np.array(["A", "C", "G", "T"], dtype=object))

    def base_codes(h):
        order = np.argsort(hashes)
        return order[np.searchsorted(hashes[order], h)]

    def groups(a1, a2):
        # Eligible pairs are non-palindromic SNVs: transitions (A/G, C/T)
        # have even code sums; transversions (A/C, G/T) have odd sums.
        return (base_codes(a1) + base_codes(a2)) % 2

    ref_codes, src_codes = codes[:ri.size], codes[ri.size:]
    ref_groups = groups(ref_a1[ri], ref_a2[ri])
    src_groups = groups(src_a1[si], src_a2[si])
    # Reuse the numeric tuple layout as (ID code, allele group, chromosome code).
    ref_keys = numeric_variant_keys(ref_codes, ref_groups, np.zeros(ri.size, dtype=np.uint64))
    _, inverse, counts = np.unique(ref_keys, return_inverse=True, return_counts=True)
    duplicates = counts[inverse] > 1
    if duplicates.any():
        warnings.warn(f"reference panel: excluded {np.count_nonzero(duplicates)} slots in "
                      f"{np.count_nonzero(counts > 1)} duplicate ID key groups",
                      RuntimeWarning, stacklevel=2)
    ri = ri[~duplicates]
    ref_keys = ref_keys[~duplicates]
    rows = np.full(reference.num_snp, -1, dtype=np.int64)
    swapped = np.zeros(reference.num_snp, dtype=bool)
    if not ri.size or not si.size:
        return rows, swapped
    src_keys = numeric_variant_keys(src_codes, src_groups, np.zeros(si.size, dtype=np.uint64))
    if chromosome_scope:
        # Reference ambiguity is assessed above, without chromosome in the key.
        ref_chr = np.concatenate([np.full(s.num_snp, i) for i, s in enumerate(reference.shards)])
        src_chr = labels.get_indexer(np.concatenate([s[4] for s in sources]))
        ref_keys["a2_hash64"] = ref_chr[ri] + 1
        src_keys["a2_hash64"] = src_chr[si] + 1
    order = np.argsort(src_keys, kind="mergesort")
    sorted_keys = src_keys[order]
    left = np.searchsorted(sorted_keys, ref_keys, side="left")
    right = np.searchsorted(sorted_keys, ref_keys, side="right")
    if np.any(right - left > 1):
        slot = ri[np.flatnonzero(right - left > 1)[0]]
        contexts = ", ".join(str(s[0]) for s in sources)
        raise ValueError(f"Ambiguous duplicate source/reference ID matching key for "
                         f"ID {ref_ids[slot]!r} in {contexts}")
    matched = right > left
    rows[ri[matched]] = si[order[left[matched]]]
    complemented = np.zeros(reference.num_snp, dtype=bool)
    slots = np.flatnonzero(rows >= 0)
    source_first = base_codes(src_a1[rows[slots]])
    target_first = base_codes(ref_a1[slots])
    target_other = base_codes(ref_a2[slots])
    swapped[slots] = (source_first == target_other) | (3 - source_first == target_other)
    complemented[slots] = (source_first != target_first) & (source_first != target_other)
    offset = 0
    for context, ids, *_ in sources:
        selected = (rows >= offset) & (rows < offset + len(ids))
        for mask, label in ((swapped, "inversion"), (complemented, "strand complementation")):
            count = np.count_nonzero(selected & mask)
            if count:
                info(f"{context}: {count} matched row(s) requiring {label}")
        offset += len(ids)
    return rows, swapped
