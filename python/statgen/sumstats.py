import json
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.errors import ParserError

from ._utils import validate_requested_shards

_CACHE_SCHEMA = "sumstats_cache/0.1"
_REQUIRED_COLS = ("chr", "bp", "a1", "a2", "z", "n")
_OPTIONAL_COLS = ("p", "beta", "se", "eaf", "info")


def _key_from_arrays(chr_arr, bp_arr, a1_arr, a2_arr) -> pd.Series:
    return (
        pd.Series(chr_arr, dtype="string")
        .str.cat(pd.Series(bp_arr, dtype=np.int64).astype("string"), sep=":")
        .str.cat(pd.Series(a1_arr, dtype="string"), sep=":")
        .str.cat(pd.Series(a2_arr, dtype="string"), sep=":")
    )


def _parse_sumstats(path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            compression="infer",
            dtype=str,
            keep_default_na=False,
            na_filter=False,
        )
    except ParserError as exc:
        raise ValueError(f"{path}: malformed TSV") from exc

    missing = [c for c in _REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing required columns: {', '.join(missing)}")

    out = df.copy()
    bp_num = pd.to_numeric(out["bp"], errors="coerce")
    bad_bp = bp_num.isna() | (np.floor(bp_num) != bp_num)
    if bad_bp.any():
        idx = int(bad_bp.idxmax())
        raise ValueError(f"{path}: row {idx + 2}: bp is not an integer: {out['bp'].iat[idx]!r}")
    out["bp"] = bp_num.astype(np.int64)

    for required_num in ("z", "n"):
        vals = pd.to_numeric(out[required_num], errors="coerce")
        bad = vals.isna() | ~np.isfinite(vals.to_numpy(dtype=float))
        if bad.any():
            idx = int(bad.idxmax())
            raise ValueError(
                f"{path}: row {idx + 2}: {required_num} must be finite numeric: {out[required_num].iat[idx]!r}"
            )
        out[required_num] = vals.astype(float)

    for col in ("chr", "a1", "a2"):
        bad_empty = out[col].eq("")
        if bad_empty.any():
            idx = int(bad_empty.idxmax())
            raise ValueError(f"{path}: row {idx + 2}: {col} must be non-empty")

    for opt in _OPTIONAL_COLS:
        if opt in out.columns:
            out[opt] = pd.to_numeric(out[opt], errors="coerce").astype(float)

    return out


class SumstatsShard:
    def __init__(
        self,
        label: str,
        checksum: str,
        zvec,
        nvec,
        logpvec,
        beta_vec=None,
        se_vec=None,
        eaf_vec=None,
        info_vec=None,
    ):
        self._label = label
        self._checksum = checksum
        self._zvec = np.asarray(zvec, dtype=float)
        self._nvec = np.asarray(nvec, dtype=float)
        self._logpvec = np.asarray(logpvec, dtype=float)
        self._beta_vec = None if beta_vec is None else np.asarray(beta_vec, dtype=float)
        self._se_vec = None if se_vec is None else np.asarray(se_vec, dtype=float)
        self._eaf_vec = None if eaf_vec is None else np.asarray(eaf_vec, dtype=float)
        self._info_vec = None if info_vec is None else np.asarray(info_vec, dtype=float)

    @property
    def label(self) -> str:
        return self._label

    @property
    def checksum(self) -> str:
        return self._checksum

    @property
    def num_snp(self) -> int:
        return self._zvec.size

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
        checksum,
        zvec,
        nvec,
        logpvec,
        beta_vec,
        se_vec,
        eaf_vec,
        info_vec,
    ):
        return cls(
            label=label,
            checksum=checksum,
            zvec=zvec,
            nvec=nvec,
            logpvec=logpvec,
            beta_vec=beta_vec,
            se_vec=se_vec,
            eaf_vec=eaf_vec,
            info_vec=info_vec,
        )


class Sumstats:
    def __init__(self, shards: list[SumstatsShard]):
        self._shards = list(shards)
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
        if not self._shards:
            return np.array([], dtype=float)
        return np.concatenate([s.zvec for s in self._shards])

    @property
    def nvec(self):
        if not self._shards:
            return np.array([], dtype=float)
        return np.concatenate([s.nvec for s in self._shards])

    @property
    def logpvec(self):
        if not self._shards:
            return np.array([], dtype=float)
        return np.concatenate([s.logpvec for s in self._shards])

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


def _derive_logp(p_vals: np.ndarray) -> np.ndarray:
    logp = np.full(p_vals.size, np.nan, dtype=float)
    finite = np.isfinite(p_vals)
    in_range = finite & (p_vals >= 0.0) & (p_vals <= 1.0)
    zero_mask = in_range & (p_vals == 0.0)
    pos_mask = in_range & (p_vals > 0.0)
    logp[zero_mask] = np.inf
    logp[pos_mask] = -np.log10(p_vals[pos_mask])
    return logp


def _build_sumstats_panel(df: pd.DataFrame, reference) -> Sumstats:
    ref_chr = np.asarray(reference.chr, dtype=object)
    ref_bp = np.asarray(reference.bp, dtype=np.int64)
    ref_a1 = np.asarray(reference.a1, dtype=object)
    ref_a2 = np.asarray(reference.a2, dtype=object)
    ref_keys = _key_from_arrays(ref_chr, ref_bp, ref_a1, ref_a2)

    src_keys = _key_from_arrays(
        df["chr"].to_numpy(dtype=object),
        df["bp"].to_numpy(dtype=np.int64),
        df["a1"].to_numpy(dtype=object),
        df["a2"].to_numpy(dtype=object),
    )
    src = df.copy()
    src["_key"] = src_keys

    key_to = {}
    for col in ("z", "n") + _OPTIONAL_COLS:
        if col in src.columns:
            key_to[col] = src.set_index("_key")[col]

    aligned_z = key_to["z"].reindex(ref_keys).to_numpy(dtype=float)
    aligned_n = key_to["n"].reindex(ref_keys).to_numpy(dtype=float)
    if "p" in key_to:
        aligned_p = key_to["p"].reindex(ref_keys).to_numpy(dtype=float)
        aligned_logp = _derive_logp(aligned_p)
    else:
        aligned_logp = np.full(ref_keys.size, np.nan, dtype=float)

    aligned_optional = {}
    for col in ("beta", "se", "eaf", "info"):
        if col in key_to:
            aligned_optional[col] = key_to[col].reindex(ref_keys).to_numpy(dtype=float)
        else:
            aligned_optional[col] = None

    shards = []
    for ref_shard, off in zip(reference.shards, reference.shard_offsets):
        start = int(off["start0"])
        stop = int(off["stop0"])
        shards.append(
            SumstatsShard(
                label=ref_shard.label,
                checksum=ref_shard.checksum,
                zvec=aligned_z[start:stop],
                nvec=aligned_n[start:stop],
                logpvec=aligned_logp[start:stop],
                beta_vec=None if aligned_optional["beta"] is None else aligned_optional["beta"][start:stop],
                se_vec=None if aligned_optional["se"] is None else aligned_optional["se"][start:stop],
                eaf_vec=None if aligned_optional["eaf"] is None else aligned_optional["eaf"][start:stop],
                info_vec=None if aligned_optional["info"] is None else aligned_optional["info"][start:stop],
            )
        )
    return Sumstats(shards)


def load_sumstats(path, reference) -> Sumstats:
    df = _parse_sumstats(Path(path))
    return _build_sumstats_panel(df, reference)


def save_sumstats_cache(sumstats: Sumstats, path) -> None:
    meta = {
        "schema": _CACHE_SCHEMA,
        "shard_labels": [s.label for s in sumstats.shards],
        "shard_checksums": [s.checksum for s in sumstats.shards],
        "has_beta": sumstats.beta_vec is not None,
        "has_se": sumstats.se_vec is not None,
        "has_eaf": sumstats.eaf_vec is not None,
        "has_info": sumstats.info_vec is not None,
    }
    arrays = {"_meta": np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8)}
    for i, s in enumerate(sumstats.shards):
        p = f"s{i}_"
        arrays[p + "zvec"] = s.zvec
        arrays[p + "nvec"] = s.nvec
        arrays[p + "logpvec"] = s.logpvec
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
                    checksum=checksums[i],
                    zvec=data[p + "zvec"],
                    nvec=data[p + "nvec"],
                    logpvec=data[p + "logpvec"],
                    beta_vec=data[p + "beta_vec"] if has_beta else None,
                    se_vec=data[p + "se_vec"] if has_se else None,
                    eaf_vec=data[p + "eaf_vec"] if has_eaf else None,
                    info_vec=data[p + "info_vec"] if has_info else None,
                )
            )
    return Sumstats(shard_objs)
