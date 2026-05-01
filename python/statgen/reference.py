import hashlib
import json
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.errors import ParserError

logger = logging.getLogger(__name__)

_CACHE_SCHEMA = "reference_cache/0.1"
_BIM_NCOLS = 6
_CANONICAL_CHR_ORDER = [str(i) for i in range(1, 23)] + ["X"]
_CHR_RANK = {c: i for i, c in enumerate(_CANONICAL_CHR_ORDER)}
_IGNORED_CHR = {"Y", "MT"}
_DNA_ALLELE_RE = re.compile(r"^[ACGT]+$")


def _validate_requested_shards(shards, available_labels, where: str) -> list[str]:
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
        if label not in _CHR_RANK:
            raise ValueError(f"{where}: unsupported shard label {label!r}; expected canonical labels 1-22 or X")
        if label in seen:
            raise ValueError(f"{where}: duplicate shard label {label!r} in shards")
        rank = _CHR_RANK[label]
        if prev_rank is not None and rank <= prev_rank:
            raise ValueError(f"{where}: shards must be in canonical subsequence order")
        seen.add(label)
        prev_rank = rank
        if label not in labels:
            raise ValueError(f"{where}: requested shard {label!r} is not present")
    return requested


def _validate_reference_sort_order(
    chr_col: pd.Series,
    bp_col: np.ndarray,
    a1_col: pd.Series,
    a2_col: pd.Series,
    line_numbers: np.ndarray,
    path: Path,
) -> None:
    chr_rank = chr_col.map(_CHR_RANK)
    bad_chr = chr_rank.isna()
    if bad_chr.any():
        idx = int(bad_chr.idxmax())
        raise ValueError(
            f"{path}:{idx + 1}: chr must use canonical labels 1-22 or X"
        )

    order_df = pd.DataFrame(
        {
            "chr_rank": chr_rank.to_numpy(dtype=np.int64),
            "bp": np.asarray(bp_col, dtype=np.int64),
            "a1": a1_col.to_numpy(dtype=object),
            "a2": a2_col.to_numpy(dtype=object),
            "line": np.asarray(line_numbers, dtype=np.int64),
        }
    )

    dup_mask = order_df.duplicated(
        subset=["chr_rank", "bp", "a1", "a2"], keep="first"
    )
    if dup_mask.any():
        line = int(order_df.loc[dup_mask.idxmax(), "line"])
        raise ValueError(
            f"{path}:{line}: duplicate (chr, bp, a1, a2) tuple is not allowed"
        )

    prev = order_df.shift(1)
    bad_order = (
        (order_df["chr_rank"] < prev["chr_rank"])
        | (
            (order_df["chr_rank"] == prev["chr_rank"])
            & (
                (order_df["bp"] < prev["bp"])
                | (
                    (order_df["bp"] == prev["bp"])
                    & (
                        (order_df["a1"] < prev["a1"])
                        | (
                            (order_df["a1"] == prev["a1"])
                            & (order_df["a2"] < prev["a2"])
                        )
                    )
                )
            )
        )
    )
    bad_order = bad_order.fillna(False)
    if bad_order.any():
        line = int(order_df.loc[bad_order.idxmax(), "line"])
        raise ValueError(
            f"{path}:{line}: rows must be sorted by (chr_rank, bp, a1, a2) in canonical contig order"
        )


def _parse_bim(path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            dtype=str,
            keep_default_na=False,
            na_filter=False,
        )
    except ParserError as exc:
        raise ValueError(
            f"{path}: expected 6 tab-separated columns"
        ) from exc

    if df.shape[1] != _BIM_NCOLS:
        raise ValueError(
            f"{path}: expected 6 tab-separated columns, got {df.shape[1]}"
        )

    chr_col = df[0]
    snp_col = df[1]
    cm_raw = df[2]
    bp_raw = df[3]
    a1_col = df[4]
    a2_col = df[5]

    bad_chr = chr_col.eq("")
    if bad_chr.any():
        idx = int(bad_chr.idxmax())
        raise ValueError(f"{path}:{idx + 1}: chr must be non-empty")

    chr_style = chr_col.str.lower().str.startswith("chr")
    if chr_style.any():
        idx = int(chr_style.idxmax())
        raise ValueError(f"{path}:{idx + 1}: chr-style labels (e.g., chr1/chrX) are not allowed")

    known_chr = chr_col.isin(_CANONICAL_CHR_ORDER) | chr_col.isin(_IGNORED_CHR)
    if not known_chr.all():
        idx = int((~known_chr).idxmax())
        raise ValueError(
            f"{path}:{idx + 1}: unsupported chr label {chr_col.iat[idx]!r}; expected 1-22, X (Y/MT are ignored)"
        )

    bad_allele = a1_col.eq("") | a2_col.eq("")
    if bad_allele.any():
        idx = int(bad_allele.idxmax())
        raise ValueError(f"{path}:{idx + 1}: a1 and a2 must be non-empty")

    bad_a1_syntax = ~a1_col.str.fullmatch(_DNA_ALLELE_RE)
    if bad_a1_syntax.any():
        idx = int(bad_a1_syntax.idxmax())
        raise ValueError(
            f"{path}:{idx + 1}: a1 must be uppercase DNA bases (A/C/G/T): {a1_col.iat[idx]!r}"
        )
    bad_a2_syntax = ~a2_col.str.fullmatch(_DNA_ALLELE_RE)
    if bad_a2_syntax.any():
        idx = int(bad_a2_syntax.idxmax())
        raise ValueError(
            f"{path}:{idx + 1}: a2 must be uppercase DNA bases (A/C/G/T): {a2_col.iat[idx]!r}"
        )

    cm_num = pd.to_numeric(cm_raw, errors="coerce")
    bad_cm = cm_num.isna()
    if bad_cm.any():
        idx = int(bad_cm.idxmax())
        raise ValueError(
            f"{path}:{idx + 1}: cm is not a number: {cm_raw.iat[idx]!r}"
        )

    bp_num = pd.to_numeric(bp_raw, errors="coerce")
    bad_bp = bp_num.isna() | (np.floor(bp_num) != bp_num)
    if bad_bp.any():
        idx = int(bad_bp.idxmax())
        raise ValueError(
            f"{path}:{idx + 1}: bp is not an integer: {bp_raw.iat[idx]!r}"
        )

    bp_int = bp_num.astype(np.int64)
    keep = chr_col.isin(_CANONICAL_CHR_ORDER)
    chr_kept = chr_col[keep]
    bp_kept = bp_int[keep]
    a1_kept = a1_col[keep]
    a2_kept = a2_col[keep]
    line_numbers = chr_kept.index.to_numpy(dtype=np.int64) + 1
    _validate_reference_sort_order(
        chr_kept,
        bp_kept.to_numpy(dtype=np.int64),
        a1_kept,
        a2_kept,
        line_numbers,
        path,
    )

    # cm is validated for BIM schema compatibility but intentionally ignored
    # in-memory per spec.
    out = pd.DataFrame(
        {
            "chr": chr_kept.to_numpy(dtype=object),
            "snp": snp_col[keep].to_numpy(dtype=object),
            "bp": bp_kept.to_numpy(dtype=np.int64),
            "a1": a1_kept.to_numpy(dtype=object),
            "a2": a2_kept.to_numpy(dtype=object),
        }
    )
    return out.reset_index(drop=True)


def _checksum_from_arrays(
    chr_arr: np.ndarray, bp_arr: np.ndarray, a1_arr: np.ndarray, a2_arr: np.ndarray
) -> str:
    if chr_arr.size == 0:
        return hashlib.md5(b"").hexdigest()

    line_series = (
        pd.Series(chr_arr, dtype="string")
        .str.cat(pd.Series(bp_arr, dtype=np.int64).astype("string"), sep=":")
        .str.cat(pd.Series(a1_arr, dtype="string"), sep=":")
        .str.cat(pd.Series(a2_arr, dtype="string"), sep=":")
    )
    text = line_series.str.cat(sep="\n") + "\n"
    return hashlib.md5(text.encode()).hexdigest()


class ReferenceShard:
    def __init__(
        self,
        label: str,
        chr_arr,
        snp_arr,
        bp_arr,
        a1_arr,
        a2_arr,
    ):
        self._label = label
        self._chr = np.asarray(chr_arr, dtype=object)
        self._snp = np.asarray(snp_arr, dtype=object)
        self._bp = np.asarray(bp_arr, dtype=np.int64)
        self._a1 = np.asarray(a1_arr, dtype=object)
        self._a2 = np.asarray(a2_arr, dtype=object)
        self._checksum = _checksum_from_arrays(self._chr, self._bp, self._a1, self._a2)

    @property
    def label(self) -> str:
        return self._label

    @property
    def num_snp(self) -> int:
        return len(self._chr)

    @property
    def chr(self) -> np.ndarray:
        return self._chr

    @property
    def snp(self) -> np.ndarray:
        return self._snp

    @property
    def bp(self) -> np.ndarray:
        return self._bp

    @property
    def a1(self) -> np.ndarray:
        return self._a1

    @property
    def a2(self) -> np.ndarray:
        return self._a2

    @property
    def checksum(self) -> str:
        return self._checksum

    @classmethod
    def _from_arrays(cls, label, chr_arr, snp_arr, bp_arr, a1_arr, a2_arr, checksum):
        obj = cls.__new__(cls)
        obj._label = label
        obj._chr = np.asarray(chr_arr, dtype=object)
        obj._snp = np.asarray(snp_arr, dtype=object)
        obj._bp = np.asarray(bp_arr, dtype=np.int64)
        obj._a1 = np.asarray(a1_arr, dtype=object)
        obj._a2 = np.asarray(a2_arr, dtype=object)
        obj._checksum = checksum
        return obj


class ReferencePanel:
    def __init__(self, shards: list):
        self._shards = list(shards)
        offsets = []
        pos = 0
        for s in self._shards:
            offsets.append({"shard_label": s.label, "start0": pos, "stop0": pos + s.num_snp})
            pos += s.num_snp
        self._shard_offsets = offsets
        self._num_snp = pos

    @property
    def num_snp(self) -> int:
        return self._num_snp

    @property
    def chr(self) -> np.ndarray:
        if not self._shards:
            return np.array([], dtype=object)
        return np.concatenate([s.chr for s in self._shards])

    @property
    def snp(self) -> np.ndarray:
        if not self._shards:
            return np.array([], dtype=object)
        return np.concatenate([s.snp for s in self._shards])

    @property
    def bp(self) -> np.ndarray:
        if not self._shards:
            return np.array([], dtype=np.int64)
        return np.concatenate([s.bp for s in self._shards])

    @property
    def a1(self) -> np.ndarray:
        if not self._shards:
            return np.array([], dtype=object)
        return np.concatenate([s.a1 for s in self._shards])

    @property
    def a2(self) -> np.ndarray:
        if not self._shards:
            return np.array([], dtype=object)
        return np.concatenate([s.a2 for s in self._shards])

    @property
    def shard_offsets(self) -> list:
        return list(self._shard_offsets)

    @property
    def shards(self) -> list:
        return list(self._shards)

    def select_shards(self, shards) -> "ReferencePanel":
        available = [s.label for s in self._shards]
        selected = _validate_requested_shards(shards, available, "ReferencePanel.select_shards")
        by_label = {s.label: s for s in self._shards}
        return ReferencePanel([by_label[label] for label in selected])

    def is_object_compatible(self, obj) -> bool:
        log_fn = logger.warning

        obj_shards = getattr(obj, "shards", None)
        if obj_shards is None:
            log_fn("statgen: is_object_compatible: object has no shards attribute")
            return False

        obj_shards = list(obj_shards)
        if len(obj_shards) != len(self._shards):
            log_fn(
                "statgen: is_object_compatible: shard count mismatch: "
                f"reference has {len(self._shards)}, object has {len(obj_shards)}"
            )
            return False

        ok = True
        for ref_s, obj_s in zip(self._shards, obj_shards):
            obj_label = getattr(obj_s, "label", None)
            if obj_label is None:
                log_fn(
                    f"statgen: is_object_compatible: shard {ref_s.label}: object shard has no label"
                )
                ok = False
                continue
            if obj_label != ref_s.label:
                log_fn(
                    "statgen: is_object_compatible: shard label mismatch: "
                    f"reference {ref_s.label}, object {obj_label}"
                )
                ok = False
                continue
            obj_num = getattr(obj_s, "num_snp", None)
            if obj_num is None:
                log_fn(
                    f"statgen: is_object_compatible: shard {ref_s.label}: "
                    "object shard has no num_snp"
                )
                ok = False
                continue
            if obj_num != ref_s.num_snp:
                log_fn(
                    f"statgen: is_object_compatible: shard {ref_s.label}: "
                    f"row count mismatch: reference {ref_s.num_snp}, object {obj_num}"
                )
                ok = False
                continue
            obj_chk = getattr(obj_s, "checksum", None)
            if obj_chk is not None and obj_chk != ref_s.checksum:
                log_fn(
                    f"statgen: is_object_compatible: shard {ref_s.label}: checksum mismatch"
                )
                ok = False
        return ok


def load_reference(path, shards=None) -> ReferencePanel:
    path_str = str(path)
    path_obj = Path(path_str)

    if "@" in path_str:
        available_labels = []
        available_paths = {}
        for label in _CANONICAL_CHR_ORDER:
            candidate = Path(path_str.replace("@", label))
            if candidate.is_file():
                available_labels.append(label)
                available_paths[label] = candidate

        if not available_labels:
            raise FileNotFoundError(
                f"No BIM shards found matching template: {path_str}"
            )

        selected = _validate_requested_shards(
            shards, available_labels, "load_reference"
        )
        shards = []
        for label in selected:
            bim_df = _parse_bim(available_paths[label])
            shards.append(
                ReferenceShard(
                    label,
                    bim_df["chr"].to_numpy(dtype=object),
                    bim_df["snp"].to_numpy(dtype=object),
                    bim_df["bp"].to_numpy(dtype=np.int64),
                    bim_df["a1"].to_numpy(dtype=object),
                    bim_df["a2"].to_numpy(dtype=object),
                )
            )
        return ReferencePanel(shards)

    bim_df = _parse_bim(path_obj)

    available_labels = [
        c for c in _CANONICAL_CHR_ORDER if (bim_df["chr"] == c).any()
    ]
    selected = _validate_requested_shards(shards, available_labels, "load_reference")
    out_shards = []
    for c in selected:
        shard_df = bim_df.loc[bim_df["chr"] == c]
        out_shards.append(
            ReferenceShard(
                c,
                shard_df["chr"].to_numpy(dtype=object),
                shard_df["snp"].to_numpy(dtype=object),
                shard_df["bp"].to_numpy(dtype=np.int64),
                shard_df["a1"].to_numpy(dtype=object),
                shard_df["a2"].to_numpy(dtype=object),
            )
        )
    return ReferencePanel(out_shards)


def save_reference_cache(panel: ReferencePanel, path) -> None:
    # Metadata (schema, labels, checksums) as a compact JSON blob stored in the npz.
    # SNP-axis numeric vectors (bp) and string arrays are stored as native
    # binary numpy arrays — not JSON — per the performance contract.
    meta = {
        "schema": _CACHE_SCHEMA,
        "shard_labels": [s.label for s in panel.shards],
        "shard_checksums": [s.checksum for s in panel.shards],
    }
    arrays: dict = {
        "_meta": np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8)
    }
    for i, s in enumerate(panel.shards):
        p = f"s{i}_"
        arrays[p + "chr"] = np.asarray(s.chr, dtype=str)
        arrays[p + "snp"] = np.asarray(s.snp, dtype=str)
        arrays[p + "bp"]  = s.bp
        arrays[p + "a1"]  = np.asarray(s.a1, dtype=str)
        arrays[p + "a2"]  = np.asarray(s.a2, dtype=str)
    np.savez_compressed(path, **arrays)


def load_reference_cache(path, shards=None) -> ReferencePanel:
    with np.load(path, allow_pickle=False) as data:
        meta = json.loads(bytes(data["_meta"]).decode())
        schema = meta.get("schema")
        if schema != _CACHE_SCHEMA:
            raise ValueError(f"Unsupported reference cache schema: {schema!r}")

        labels = list(meta.get("shard_labels", []))
        checksums = list(meta.get("shard_checksums", []))
        if len(labels) != len(checksums):
            raise ValueError("Invalid reference cache: shard_labels and shard_checksums length mismatch")
        selected = _validate_requested_shards(shards, labels, "load_reference_cache")
        label_to_index = {label: i for i, label in enumerate(labels)}

        shard_objs = []
        for label in selected:
            i = label_to_index[label]
            p = f"s{i}_"
            shard_objs.append(
                ReferenceShard._from_arrays(
                    label,
                    data[p + "chr"],
                    data[p + "snp"],
                    data[p + "bp"],
                    data[p + "a1"],
                    data[p + "a2"],
                    checksums[i],
                )
            )
    return ReferencePanel(shard_objs)
