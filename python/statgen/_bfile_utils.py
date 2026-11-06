from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError, ParserError, ParserWarning

from ._utils import CANONICAL_CHR_ORDER, CHR_RANK, allele_hash64
from ._variant_match import check_unique_sorted_variant_keys, numeric_variant_keys

IGNORED_CHR = {"Y", "MT"}
FAM_PUBLIC_COLUMNS = ["fid", "iid", "father_id", "mother_id", "sex"]

_BIM_NAMES = ["chr", "snp", "cm", "bp", "a1", "a2"]
_FAM_NAMES = ["fid", "iid", "father_id", "mother_id", "sex", "pheno"]
_PLOIDY_NAMES = ["ploidy_male", "ploidy_female"]
_DNA_ALLELE_RE = re.compile(r"^[ACGT]+$")
_BED_LOOKUP_INT8 = None
_BED_MAX_READ_BYTES = 64 * 1024 * 1024


def _read_csv_strict(path: Path, *, sep, names, description: str | None = None) -> pd.DataFrame:
    expected = description or f"{len(names)} delimited columns"
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", ParserWarning)
            df = pd.read_csv(
                path,
                sep=sep,
                header=None,
                dtype=str,
                index_col=False,
                keep_default_na=False,
                na_filter=False,
            )
    except (EmptyDataError, ParserError, ParserWarning, ValueError) as exc:
        raise ValueError(f"{path}: expected {expected}") from exc
    if df.shape[1] != len(names):
        raise ValueError(f"{path}: expected {expected}, got {df.shape[1]}")
    df.columns = names
    return df


def _validate_source_sort_order(source_df: pd.DataFrame, path: Path) -> None:
    if source_df.empty:
        return
    line = (
        source_df["line"].to_numpy(dtype=np.int64)
        if "line" in source_df.columns
        else source_df.index.to_numpy(dtype=np.int64) + 1
    )
    chr_rank = source_df["chr"].map(CHR_RANK)
    bad_chr = chr_rank.isna()
    if bad_chr.any():
        idx = int(bad_chr.idxmax())
        raise ValueError(f"{path}:{idx + 1}: chr must use canonical labels 1-22 or X")

    keys = numeric_variant_keys(
        source_df["bp"].to_numpy(dtype=np.int64),
        allele_hash64(source_df["a1"].to_numpy(dtype=object)),
        allele_hash64(source_df["a2"].to_numpy(dtype=object)),
    )
    order_df = pd.DataFrame(
        {
            "chr_rank": chr_rank.to_numpy(dtype=np.int64),
            "bp": source_df["bp"].to_numpy(dtype=np.int64),
            "line": line,
        }
    )

    for rank in np.unique(order_df["chr_rank"].to_numpy(dtype=np.int64)):
        mask = order_df["chr_rank"].to_numpy(dtype=np.int64) == rank
        shard_keys = keys[mask]
        sorted_keys = np.sort(shard_keys, order=("bp", "a1_hash64", "a2_hash64"))
        check_unique_sorted_variant_keys(sorted_keys, f"BIM {path}", "variant")

    prev = order_df.shift(1)
    bad_order = (
        (order_df["chr_rank"] < prev["chr_rank"])
        | (
            (order_df["chr_rank"] == prev["chr_rank"])
            & (order_df["bp"] < prev["bp"])
        )
    ).fillna(False)
    if bad_order.any():
        line = int(order_df.loc[bad_order.idxmax(), "line"])
        raise ValueError(
            f"{path}:{line}: rows must be sorted by (chr_rank, bp) in canonical contig order"
        )


def _parse_bim(path: Path) -> pd.DataFrame:
    path = Path(path)
    df = _read_csv_strict(
        path,
        sep=r"\s+",
        names=_BIM_NAMES,
        description="6 whitespace-delimited columns",
    )
    df = df.copy()

    chr_col = df["chr"]
    bad_chr = chr_col.eq("")
    if bad_chr.any():
        idx = int(bad_chr.idxmax())
        raise ValueError(f"{path}:{idx + 1}: chr must be non-empty")

    chr_style = chr_col.str.lower().str.startswith("chr")
    if chr_style.any():
        idx = int(chr_style.idxmax())
        raise ValueError(f"{path}:{idx + 1}: chr-style labels (e.g., chr1/chrX) are not allowed")

    known_chr = chr_col.isin(CANONICAL_CHR_ORDER) | chr_col.isin(IGNORED_CHR)
    if not known_chr.all():
        idx = int((~known_chr).idxmax())
        raise ValueError(
            f"{path}:{idx + 1}: unsupported chr label {chr_col.iat[idx]!r}; expected 1-22, X (Y/MT are ignored)"
        )

    bad_allele = df["a1"].eq("") | df["a2"].eq("")
    if bad_allele.any():
        idx = int(bad_allele.idxmax())
        raise ValueError(f"{path}:{idx + 1}: a1 and a2 must be non-empty")

    bad_a1_syntax = ~df["a1"].str.fullmatch(_DNA_ALLELE_RE.pattern)
    if bad_a1_syntax.any():
        idx = int(bad_a1_syntax.idxmax())
        raise ValueError(
            f"{path}:{idx + 1}: a1 must be uppercase DNA bases (A/C/G/T): {df['a1'].iat[idx]!r}"
        )
    bad_a2_syntax = ~df["a2"].str.fullmatch(_DNA_ALLELE_RE.pattern)
    if bad_a2_syntax.any():
        idx = int(bad_a2_syntax.idxmax())
        raise ValueError(
            f"{path}:{idx + 1}: a2 must be uppercase DNA bases (A/C/G/T): {df['a2'].iat[idx]!r}"
        )
    same_allele = df["a1"].eq(df["a2"])
    if same_allele.any():
        idx = int(same_allele.idxmax())
        raise ValueError(f"{path}:{idx + 1}: a1 and a2 must differ")

    cm_num = pd.to_numeric(df["cm"], errors="coerce")
    bad_cm = cm_num.isna()
    if bad_cm.any():
        idx = int(bad_cm.idxmax())
        raise ValueError(f"{path}:{idx + 1}: cm is not a number: {df['cm'].iat[idx]!r}")

    bp_num = pd.to_numeric(df["bp"], errors="coerce")
    bad_bp = bp_num.isna() | (np.floor(bp_num) != bp_num)
    if bad_bp.any():
        idx = int(bad_bp.idxmax())
        raise ValueError(f"{path}:{idx + 1}: bp is not an integer: {df['bp'].iat[idx]!r}")
    df["bp"] = bp_num.astype(np.int64)

    _validate_source_sort_order(df.loc[df["chr"].isin(CANONICAL_CHR_ORDER)], path)
    return df.reset_index(drop=True)


def _parse_fam(path: Path) -> pd.DataFrame:
    path = Path(path)
    df = _read_csv_strict(
        path,
        sep=r"\s+",
        names=_FAM_NAMES,
        description="6 whitespace-delimited columns",
    )
    required = df[_FAM_NAMES].eq("")
    if required.any().any():
        row = int(np.flatnonzero(required.any(axis=1).to_numpy())[0] + 1)
        raise ValueError(f"{path}:{row}: FAM must contain exactly 6 non-empty columns")

    sex_num = pd.to_numeric(df["sex"], errors="coerce")
    bad_sex = sex_num.isna() | (np.floor(sex_num) != sex_num) | ~sex_num.isin([0, 1, 2])
    if bad_sex.any():
        idx = int(bad_sex.idxmax())
        raise ValueError(f"{path}:{idx + 1}: FAM sex must be one of PLINK values 1, 2, or 0")

    out = df[["fid", "iid", "father_id", "mother_id"]].copy()
    out["sex"] = sex_num.astype(np.int8)
    dup = out.duplicated(subset=["fid", "iid"], keep="first")
    if dup.any():
        idx = int(dup.idxmax())
        raise ValueError(f"{path}:{idx + 1}: duplicate FAM subject pair (fid, iid)")
    return out.reset_index(drop=True)


def _parse_ploidy(path: Path | None, source_num_snp: int, source_chr=None) -> tuple[np.ndarray, np.ndarray]:
    if path is None:
        ploidy_male = np.full(source_num_snp, 2.0, dtype=float)
        ploidy_female = np.full(source_num_snp, 2.0, dtype=float)
        if source_chr is not None:
            chrs = np.asarray(source_chr, dtype=object).reshape(-1)
            if chrs.size != source_num_snp:
                raise ValueError("source_chr length must match source_num_snp")
            ploidy_male[chrs == "X"] = 1.0
        return ploidy_male, ploidy_female

    path = Path(path)
    df = _read_csv_strict(
        path,
        sep="\t",
        names=_PLOIDY_NAMES,
        description="2 tab-delimited columns",
    )
    if df.shape[0] != source_num_snp:
        raise ValueError(
            f"{path}: PLOIDY row count mismatch: expected {source_num_snp}, got {df.shape[0]}"
        )
    numeric = df[_PLOIDY_NAMES].apply(pd.to_numeric, errors="coerce")
    bad_numeric = numeric.isna()
    if bad_numeric.any().any():
        row = int(np.flatnonzero(bad_numeric.any(axis=1).to_numpy())[0] + 1)
        raise ValueError(f"{path}:{row}: ploidy values must be integers 0, 1, or 2")
    values = numeric.to_numpy(dtype=np.int64)
    bad = ~np.isin(values, [0, 1, 2])
    if bad.any():
        row = int(np.flatnonzero(bad.any(axis=1))[0] + 1)
        raise ValueError(f"{path}:{row}: ploidy values must be 0, 1, or 2")
    return (
        numeric["ploidy_male"].to_numpy(dtype=float),
        numeric["ploidy_female"].to_numpy(dtype=float),
    )


def _expected_bed_size(source_num_sample: int, source_num_snp: int) -> int:
    return 3 + ((int(source_num_sample) + 3) // 4) * int(source_num_snp)


def _validate_bed(path: Path, source_num_sample: int, source_num_snp: int) -> int:
    path = Path(path)
    with open(path, "rb") as f:
        magic = f.read(3)
    if len(magic) != 3 or magic[:2] != b"\x6c\x1b":
        raise ValueError(f"{path}: invalid PLINK BED magic bytes")
    if magic[2] != 1:
        raise ValueError(f"{path}: PLINK BED must be SNP-major mode")
    size = path.stat().st_size
    expected = _expected_bed_size(source_num_sample, source_num_snp)
    if size != expected:
        raise ValueError(f"{path}: BED file size mismatch: expected {expected}, got {size}")
    return int(size)


def _bed_lookup_int8() -> np.ndarray:
    global _BED_LOOKUP_INT8
    if _BED_LOOKUP_INT8 is None:
        values = np.arange(256, dtype=np.uint16)[:, None]
        shifts = (np.arange(4, dtype=np.uint16) * 2)[None, :]
        two_bit = ((values >> shifts) & 0b11).astype(np.uint8)
        lookup = np.empty((256, 4), dtype=np.int8)
        lookup[two_bit == 0b00] = 2
        lookup[two_bit == 0b01] = -1
        lookup[two_bit == 0b10] = 1
        lookup[two_bit == 0b11] = 0
        lookup.setflags(write=False)
        _BED_LOOKUP_INT8 = lookup
    return _BED_LOOKUP_INT8


def _read_bed_rows_int8(path: Path, source_rows: np.ndarray, source_num_sample: int) -> np.ndarray:
    path = Path(path)
    source_rows = np.asarray(source_rows, dtype=np.int64).reshape(-1)
    n_samples = int(source_num_sample)
    bytes_per_snp = (n_samples + 3) // 4
    if source_rows.size == 0:
        return np.empty((n_samples, 0), dtype=np.int8)
    unique_rows, inverse = np.unique(source_rows, return_inverse=True)
    unique_decoded = np.empty((n_samples, unique_rows.size), dtype=np.int8)
    rows_per_read = unique_rows.size if bytes_per_snp == 0 else max(1, _BED_MAX_READ_BYTES // bytes_per_snp)
    lookup = _bed_lookup_int8()
    with open(path, "rb") as f:
        group_starts = np.r_[0, np.flatnonzero(np.diff(unique_rows) != 1) + 1]
        group_stops = np.r_[group_starts[1:], unique_rows.size]
        # PERF: one seek/read per contiguous BED row block; requested order and
        #       repeats are restored after block decoding. Blocks are capped to
        #       bound transient packed/decode allocations for large requests.
        for start, stop in zip(group_starts, group_stops):
            for chunk_start in range(start, stop, rows_per_read):
                chunk_stop = min(chunk_start + rows_per_read, stop)
                block_rows = unique_rows[chunk_start:chunk_stop]
                expected_bytes = int(bytes_per_snp * block_rows.size)
                f.seek(3 + int(block_rows[0]) * bytes_per_snp)
                packed = f.read(expected_bytes)
                if len(packed) != expected_bytes:
                    raise ValueError(
                        f"{path}: short BED read for source SNP rows {int(block_rows[0])}-{int(block_rows[-1])}: "
                        f"expected {expected_bytes} bytes, got {len(packed)}"
                    )
                byte_values = np.frombuffer(packed, dtype=np.uint8)
                block_decoded = lookup[byte_values].reshape(block_rows.size, bytes_per_snp * 4)
                unique_decoded[:, chunk_start:chunk_stop] = block_decoded[:, :n_samples].T
    return unique_decoded[:, inverse]
