import json
import re
import warnings
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError, ParserError

from ._utils import CANONICAL_CHR_ORDER, _validate_allow_indel_swaps, allele_hash64, validate_requested_shards
from ._verbosity import report_variant_swaps
from ._variant_match import match_shard_numeric, numeric_variant_keys

_CACHE_SCHEMA = "sumstats_cache/0.1"
_IGNORED_CHR = {"Y", "MT"}
_DNA_ALLELE_RE = re.compile(r"^[ACGT]+$")
_REQUIRED_COLS = ("chr", "bp", "a1", "a2", "p")
_OPTIONAL_COLS = ("z", "n", "beta", "se", "eaf", "info")
_OPTIONAL_VECTOR_FIELDS = ("z", "n", "beta", "se", "eaf", "info")
_SUMSTATS_COL_MAP = {
    "pos": "bp",
    "effectallele": "a1",
    "otherallele": "a2",
}
_NUMERIC_NA_VALUES = ["", "NaN", "nan", "NA", "N/A", "na", "NULL", "null"]


def _canonicalize_column_names(columns, path: Path) -> list[str]:
    columns = [_SUMSTATS_COL_MAP.get(str(c).lower(), str(c).lower()) for c in columns]
    duplicates = sorted({c for c in columns if columns.count(c) > 1})
    if duplicates:
        raise ValueError(f"{path}: duplicate columns after column normalization: {', '.join(duplicates)}")
    return columns


def _canonicalize_columns(df: pd.DataFrame, path: Path) -> pd.DataFrame:
    out = df.copy()
    out.columns = _canonicalize_column_names(df.columns, path)
    return out


def _read_sumstats_header(path: Path) -> list[str]:
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    try:
        with opener(path, "rt", newline="") as f:
            header = f.readline()
    except OSError as exc:
        raise ValueError(f"{path}: malformed TSV") from exc
    if header == "":
        raise ValueError(f"{path}: malformed TSV")
    return header.rstrip("\r\n").split("\t")


def _parse_sumstats(path: Path) -> pd.DataFrame:
    raw_columns = _read_sumstats_header(path)
    columns = _canonicalize_column_names(raw_columns, path)
    missing = [c for c in _REQUIRED_COLS if c not in columns]
    if missing:
        raise ValueError(f"{path}: missing required columns: {', '.join(missing)}")

    selected = [c in (*_REQUIRED_COLS, *_OPTIONAL_COLS) for c in columns]
    selected_raw = [c for c, keep in zip(raw_columns, selected) if keep]
    selected_names = [c for c, keep in zip(columns, selected) if keep]
    char_cols = {
        raw: str
        for raw, canonical in zip(selected_raw, selected_names)
        if canonical in {"chr", "a1", "a2"}
    }
    numeric_cols = {
        raw: float
        for raw, canonical in zip(selected_raw, selected_names)
        if canonical not in {"chr", "a1", "a2"}
    }
    numeric_na_values = {raw: _NUMERIC_NA_VALUES for raw in numeric_cols}

    try:
        df = pd.read_csv(
            path,
            sep="\t",
            compression="infer",
            usecols=selected_raw,
            dtype={**char_cols, **numeric_cols},
            keep_default_na=False,
            na_values=numeric_na_values,
        )
    except (EmptyDataError, ParserError, ValueError) as exc:
        raise ValueError(f"{path}: malformed TSV") from exc

    df.columns = selected_names

    missing = [c for c in _REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing required columns: {', '.join(missing)}")

    out = df
    for col in ("chr", "a1", "a2"):
        bad_empty = out[col].eq("")
        if bad_empty.any():
            idx = int(bad_empty.idxmax())
            raise ValueError(f"{path}: row {idx + 2}: {col} must be non-empty")

    chr_col = out["chr"]
    chr_style = chr_col.str.lower().str.startswith("chr")
    if chr_style.any():
        idx = int(chr_style.idxmax())
        raise ValueError(f"{path}: row {idx + 2}: chr-style labels (e.g., chr1/chrX) are not allowed")
    known_chr = chr_col.isin(CANONICAL_CHR_ORDER) | chr_col.isin(_IGNORED_CHR)
    if not known_chr.all():
        idx = int((~known_chr).idxmax())
        raise ValueError(
            f"{path}: row {idx + 2}: unsupported chr label {chr_col.iat[idx]!r}; expected 1-22, X (Y/MT are ignored)"
        )
    out = out.loc[chr_col.isin(CANONICAL_CHR_ORDER)].copy()

    bad_a1_syntax = ~out["a1"].str.fullmatch(_DNA_ALLELE_RE.pattern)
    if bad_a1_syntax.any():
        idx = int(bad_a1_syntax.idxmax())
        raise ValueError(
            f"{path}: row {idx + 2}: a1 must be uppercase DNA bases (A/C/G/T): {out.loc[idx, 'a1']!r}"
        )
    bad_a2_syntax = ~out["a2"].str.fullmatch(_DNA_ALLELE_RE.pattern)
    if bad_a2_syntax.any():
        idx = int(bad_a2_syntax.idxmax())
        raise ValueError(
            f"{path}: row {idx + 2}: a2 must be uppercase DNA bases (A/C/G/T): {out.loc[idx, 'a2']!r}"
        )
    same_allele = out["a1"].eq(out["a2"])
    if same_allele.any():
        idx = int(same_allele.idxmax())
        raise ValueError(f"{path}: row {idx + 2}: a1 and a2 must differ")

    bad_bp = out["bp"].isna() | (np.floor(out["bp"]) != out["bp"])
    if bad_bp.any():
        idx = int(bad_bp.idxmax())
        raise ValueError(f"{path}: row {idx + 2}: bp is not an integer: {out.loc[idx, 'bp']!r}")
    out["bp"] = out["bp"].astype(np.int64)

    p_arr = out["p"].to_numpy(dtype=float)
    bad_p = out["p"].isna() | ~np.isfinite(p_arr) | (p_arr < 0.0) | (p_arr > 1.0)
    if bad_p.any():
        idx = int(bad_p.idxmax())
        raise ValueError(f"{path}: row {idx + 2}: p must be finite numeric in [0, 1]: {out.loc[idx, 'p']!r}")

    for opt in _OPTIONAL_COLS:
        if opt in out.columns:
            out[opt] = out[opt].mask(~np.isfinite(out[opt].to_numpy(dtype=float)))

    return out


class SumstatsShard:
    def __init__(
        self,
        label: str,
        reference_checksum: str,
        logpvec,
        zvec=None,
        nvec=None,
        beta_vec=None,
        se_vec=None,
        eaf_vec=None,
        info_vec=None,
    ):
        self._label = label
        self._reference_checksum = str(reference_checksum)
        self._logpvec = np.asarray(logpvec, dtype=float).reshape(-1)
        n = self._logpvec.size
        self._zvec = None if zvec is None else _as_optional_shard_vec("zvec", zvec, n)
        self._nvec = None if nvec is None else _as_optional_shard_vec("nvec", nvec, n)
        self._beta_vec = None if beta_vec is None else np.asarray(beta_vec, dtype=float).reshape(-1)
        self._se_vec = None if se_vec is None else np.asarray(se_vec, dtype=float).reshape(-1)
        self._eaf_vec = None if eaf_vec is None else np.asarray(eaf_vec, dtype=float).reshape(-1)
        self._info_vec = None if info_vec is None else np.asarray(info_vec, dtype=float).reshape(-1)
        for name in ("beta_vec", "se_vec", "eaf_vec", "info_vec"):
            vec = getattr(self, f"_{name}")
            if vec is not None and vec.reshape(-1).size != n:
                raise ValueError(f"{name} length mismatch: expected {n}, got {vec.reshape(-1).size}")

    @property
    def label(self) -> str:
        return self._label

    @property
    def reference_checksum(self) -> str:
        return self._reference_checksum

    @property
    def num_snp(self) -> int:
        return self._logpvec.size

    @property
    def zvec(self) -> np.ndarray:
        return self._zvec

    @property
    def nvec(self) -> np.ndarray:
        return self._nvec

    @property
    def logpvec(self) -> np.ndarray:
        return self._logpvec

    @property
    def is_present(self) -> np.ndarray:
        return ~np.isnan(self._logpvec)

    @property
    def beta_vec(self):
        return self._beta_vec

    @property
    def se_vec(self):
        return self._se_vec

    @property
    def eaf_vec(self):
        return self._eaf_vec

    @property
    def info_vec(self):
        return self._info_vec

    @classmethod
    def _from_arrays(
        cls,
        label,
        reference_checksum,
        logpvec,
        zvec,
        nvec,
        beta_vec,
        se_vec,
        eaf_vec,
        info_vec,
    ):
        return cls(
            label=label,
            reference_checksum=reference_checksum,
            logpvec=logpvec,
            zvec=zvec,
            nvec=nvec,
            beta_vec=beta_vec,
            se_vec=se_vec,
            eaf_vec=eaf_vec,
            info_vec=info_vec,
        )


class Sumstats:
    def __init__(self, shards: list[SumstatsShard]):
        self._shards = list(shards)
        if not self._shards:
            raise ValueError("Sumstats requires at least one shard")
        self._num_snp = int(sum(s.num_snp for s in self._shards))
        self._shard_offsets = []
        pos = 0
        for shard in self._shards:
            self._shard_offsets.append(
                {"shard_label": shard.label, "start0": pos, "stop0": pos + shard.num_snp}
            )
            pos += shard.num_snp

    @property
    def shards(self) -> list[SumstatsShard]:
        return list(self._shards)

    @property
    def num_snp(self) -> int:
        return self._num_snp

    @property
    def shard_offsets(self) -> list:
        return list(self._shard_offsets)

    @property
    def zvec(self):
        return self._optional_concat("zvec")

    @property
    def nvec(self):
        return self._optional_concat("nvec")

    @property
    def logpvec(self):
        if not self._shards:
            return np.array([], dtype=float)
        return np.concatenate([s.logpvec for s in self._shards])

    @property
    def is_present(self):
        return ~np.isnan(self.logpvec)

    def _optional_concat(self, field_name: str):
        if not self._shards:
            return None
        first = getattr(self._shards[0], field_name)
        if first is None:
            return None
        return np.concatenate([getattr(s, field_name) for s in self._shards])

    @property
    def beta_vec(self):
        return self._optional_concat("beta_vec")

    @property
    def se_vec(self):
        return self._optional_concat("se_vec")

    @property
    def eaf_vec(self):
        return self._optional_concat("eaf_vec")

    @property
    def info_vec(self):
        return self._optional_concat("info_vec")

    def select_shards(self, shards) -> "Sumstats":
        available = [s.label for s in self._shards]
        selected = validate_requested_shards(shards, available, "Sumstats.select_shards")
        by_label = {s.label: s for s in self._shards}
        return Sumstats([by_label[label] for label in selected])

    def save_cache(self, path) -> None:
        save_sumstats_cache(self, path)


def _derive_logp(p_vals: np.ndarray) -> np.ndarray:
    logp = np.full(p_vals.size, np.nan, dtype=float)
    finite = np.isfinite(p_vals)
    in_range = finite & (p_vals >= 0.0) & (p_vals <= 1.0)
    zero_mask = in_range & (p_vals == 0.0)
    pos_mask = in_range & (p_vals > 0.0)
    logp[zero_mask] = np.inf
    logp[pos_mask] = -np.log10(p_vals[pos_mask])
    return logp


def _validate_pvec(p_vals: np.ndarray, name: str) -> None:
    bad = ~np.isfinite(p_vals) | (p_vals < 0.0) | (p_vals > 1.0)
    if bad.any():
        idx = int(np.flatnonzero(bad)[0])
        raise ValueError(f"{name}[{idx}] must be finite numeric in [0, 1]")


def _as_optional_shard_vec(name: str, vec, n: int) -> np.ndarray:
    arr = np.asarray(vec, dtype=float).reshape(-1)
    if arr.size != n:
        raise ValueError(f"{name} length mismatch: expected {n}, got {arr.size}")
    return arr


def _coerce_aligned_vec(name: str, vec, n: int, allow_inf: bool) -> np.ndarray:
    arr = np.asarray(vec, dtype=float)
    if arr.ndim > 2 or (arr.ndim == 2 and 1 not in arr.shape):
        raise ValueError(f"{name} must be a vector with length {n}")
    arr = arr.reshape(-1)
    if arr.size != n:
        raise ValueError(f"{name} length mismatch: expected {n}, got {arr.size}")
    if not allow_inf:
        bad = ~(np.isfinite(arr) | np.isnan(arr))
        if bad.any():
            idx = int(np.flatnonzero(bad)[0])
            raise ValueError(f"{name}[{idx}] must be finite numeric or NaN")
    return arr


def _build_sumstats_from_aligned(
    reference,
    aligned_logp: np.ndarray,
    aligned_optional: dict[str, np.ndarray | None],
) -> Sumstats:
    n = int(reference.num_snp)
    if aligned_logp.size != n:
        raise ValueError("aligned vector length mismatch with reference")
    for field, vec in aligned_optional.items():
        if vec is not None and vec.size != n:
            raise ValueError(f"aligned {field} vector length mismatch with reference")

    shards = []
    for ref_shard, off in zip(reference.shards, reference.shard_offsets):
        start = int(off["start0"])
        stop = int(off["stop0"])
        shards.append(
            SumstatsShard(
                label=ref_shard.label,
                reference_checksum=ref_shard.checksum,
                logpvec=aligned_logp[start:stop],
                zvec=None if aligned_optional["z"] is None else aligned_optional["z"][start:stop],
                nvec=None if aligned_optional["n"] is None else aligned_optional["n"][start:stop],
                beta_vec=None if aligned_optional["beta"] is None else aligned_optional["beta"][start:stop],
                se_vec=None if aligned_optional["se"] is None else aligned_optional["se"][start:stop],
                eaf_vec=None if aligned_optional["eaf"] is None else aligned_optional["eaf"][start:stop],
                info_vec=None if aligned_optional["info"] is None else aligned_optional["info"][start:stop],
            )
        )
    return Sumstats(shards)


def _build_sumstats_panel(df: pd.DataFrame, reference, source_path, allow_indel_swaps=False) -> Sumstats:
    n = int(reference.num_snp)
    aligned_optional = {
        col: (np.full(n, np.nan, dtype=float) if col in df.columns else None)
        for col in _OPTIONAL_VECTOR_FIELDS
    }
    aligned_p = np.full(n, np.nan, dtype=float)

    src_chr = df["chr"].to_numpy(dtype=object)
    src_bp = df["bp"].to_numpy(dtype=np.int64)
    src_a1_hash64 = allele_hash64(df["a1"].to_numpy(dtype=object))
    src_a2_hash64 = allele_hash64(df["a2"].to_numpy(dtype=object))

    for ref_shard, off in zip(reference.shards, reference.shard_offsets):
        start = int(off["start0"])
        stop = int(off["stop0"])
        src_mask = src_chr == ref_shard.label
        src_idx = np.flatnonzero(src_mask)
        ref_keys = numeric_variant_keys(ref_shard.bp, ref_shard.a1_hash64, ref_shard.a2_hash64)
        src_keys = numeric_variant_keys(src_bp[src_idx], src_a1_hash64[src_idx], src_a2_hash64[src_idx])
        match = match_shard_numeric(
            ref_keys, src_keys, ref_shard.label, source_name="sumstats",
            allow_indel_swaps=allow_indel_swaps,
        )
        local_match = match.source_rows
        report_variant_swaps(source_path, ref_shard.label, np.count_nonzero(match.is_allele_swapped),
                             match.num_strand_ambiguous, match.num_disabled_indels)
        has_match = local_match >= 0
        matched_src = src_idx[local_match[has_match]]
        aligned_ix = np.arange(start, stop, dtype=np.int64)[has_match]

        aligned_p[aligned_ix] = df["p"].to_numpy(dtype=float)[matched_src]
        for col, out in aligned_optional.items():
            if out is not None:
                out[aligned_ix] = df[col].to_numpy(dtype=float)[matched_src]
        swapped_ix = np.arange(start, stop)[match.is_allele_swapped]
        for col in ("z", "beta", "eaf"):
            out = aligned_optional[col]
            if out is not None:
                out[swapped_ix] = 1 - out[swapped_ix] if col == "eaf" else -out[swapped_ix]
        if not has_match.any():
            warnings.warn(f"{source_path}: no matching source rows for reference shard {ref_shard.label!r}",
                          RuntimeWarning, stacklevel=3)

    aligned_logp = _derive_logp(aligned_p)
    _warn_optional_zn_completeness(
        aligned_optional["z"],
        aligned_optional["n"],
        aligned_logp,
        context="load_sumstats",
        stacklevel=3,
    )

    return _build_sumstats_from_aligned(
        reference,
        aligned_logp=aligned_logp,
        aligned_optional=aligned_optional,
    )


def load_sumstats(path, reference, allow_indel_swaps=False) -> Sumstats:
    allow_indel_swaps = _validate_allow_indel_swaps(allow_indel_swaps)
    path = Path(path)
    return _build_sumstats_panel(_parse_sumstats(path), reference, path, allow_indel_swaps)


def create_sumstats(
    reference,
    p,
    z=None,
    n=None,
    beta=None,
    se=None,
    eaf=None,
    info=None,
) -> Sumstats:
    num_snp_value = int(reference.num_snp)
    aligned_p = _coerce_aligned_vec("p", p, n=num_snp_value, allow_inf=True)
    _validate_pvec(aligned_p, "p")
    aligned_logp = _derive_logp(aligned_p)

    aligned_optional = {
        "z": None if z is None else _coerce_aligned_vec("z", z, n=num_snp_value, allow_inf=False),
        "n": None if n is None else _coerce_aligned_vec("n", n, n=num_snp_value, allow_inf=False),
        "beta": None if beta is None else _coerce_aligned_vec("beta", beta, n=num_snp_value, allow_inf=False),
        "se": None if se is None else _coerce_aligned_vec("se", se, n=num_snp_value, allow_inf=False),
        "eaf": None if eaf is None else _coerce_aligned_vec("eaf", eaf, n=num_snp_value, allow_inf=False),
        "info": None if info is None else _coerce_aligned_vec("info", info, n=num_snp_value, allow_inf=False),
    }

    return _build_sumstats_from_aligned(
        reference,
        aligned_logp=aligned_logp,
        aligned_optional=aligned_optional,
    )


def save_sumstats_cache(sumstats: Sumstats, path) -> None:
    meta = {
        "schema": _CACHE_SCHEMA,
        "shard_labels": [s.label for s in sumstats.shards],
        "shard_checksums": [s.reference_checksum for s in sumstats.shards],
        "has_z": sumstats.zvec is not None,
        "has_n": sumstats.nvec is not None,
        "has_beta": sumstats.beta_vec is not None,
        "has_se": sumstats.se_vec is not None,
        "has_eaf": sumstats.eaf_vec is not None,
        "has_info": sumstats.info_vec is not None,
    }
    arrays = {"_meta": np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8)}
    for i, s in enumerate(sumstats.shards):
        p = f"s{i}_"
        arrays[p + "logpvec"] = s.logpvec
        if meta["has_z"]:
            arrays[p + "zvec"] = s.zvec
        if meta["has_n"]:
            arrays[p + "nvec"] = s.nvec
        if meta["has_beta"]:
            arrays[p + "beta_vec"] = s.beta_vec
        if meta["has_se"]:
            arrays[p + "se_vec"] = s.se_vec
        if meta["has_eaf"]:
            arrays[p + "eaf_vec"] = s.eaf_vec
        if meta["has_info"]:
            arrays[p + "info_vec"] = s.info_vec
    np.savez_compressed(path, **arrays)


def load_sumstats_cache(path, shards=None) -> Sumstats:
    with np.load(path, allow_pickle=False) as data:
        meta = json.loads(bytes(data["_meta"]).decode())
        schema = meta.get("schema")
        if schema != _CACHE_SCHEMA:
            raise ValueError(f"Unsupported sumstats cache schema: {schema!r}")

        labels = list(meta.get("shard_labels", []))
        checksums = list(meta.get("shard_checksums", []))
        if len(labels) != len(checksums):
            raise ValueError("Invalid sumstats cache: shard_labels and shard_checksums length mismatch")

        selected = validate_requested_shards(shards, labels, "load_sumstats_cache")
        label_to_index = {label: i for i, label in enumerate(labels)}

        if "has_z" not in meta or "has_n" not in meta:
            raise ValueError("Invalid sumstats cache: missing has_z/has_n metadata")
        has_z = bool(meta["has_z"])
        has_n = bool(meta["has_n"])
        has_beta = bool(meta.get("has_beta", False))
        has_se = bool(meta.get("has_se", False))
        has_eaf = bool(meta.get("has_eaf", False))
        has_info = bool(meta.get("has_info", False))

        shard_objs = []
        for label in selected:
            i = label_to_index[label]
            p = f"s{i}_"
            shard_objs.append(
                SumstatsShard._from_arrays(
                    label=label,
                    reference_checksum=checksums[i],
                    logpvec=data[p + "logpvec"],
                    zvec=data[p + "zvec"] if has_z else None,
                    nvec=data[p + "nvec"] if has_n else None,
                    beta_vec=data[p + "beta_vec"] if has_beta else None,
                    se_vec=data[p + "se_vec"] if has_se else None,
                    eaf_vec=data[p + "eaf_vec"] if has_eaf else None,
                    info_vec=data[p + "info_vec"] if has_info else None,
                )
            )
    out = Sumstats(shard_objs)
    _warn_optional_zn_completeness(
        out.zvec,
        out.nvec,
        out.logpvec,
        context="load_sumstats_cache",
        stacklevel=2,
    )
    return out


def _warn_optional_zn_completeness(zvec, nvec, logpvec, context: str, stacklevel: int) -> None:
    is_present = ~np.isnan(logpvec)
    for name, vec in (("zvec", zvec), ("nvec", nvec)):
        if vec is None:
            warnings.warn(f"{context}: {name} is absent", RuntimeWarning, stacklevel=stacklevel)
            continue
        if np.isnan(vec[is_present]).any():
            warnings.warn(
                f"{context}: {name} has missing values among present sumstats variants",
                RuntimeWarning,
                stacklevel=stacklevel,
            )
