import hashlib
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.errors import ParserError

logger = logging.getLogger(__name__)

_CACHE_SCHEMA = "reference_cache/0.1"
_BIM_NCOLS = 6
_CANONICAL_CHR_ORDER = [str(i) for i in range(1, 23)] + ["X"]
_CHR_RANK = {c: i for i, c in enumerate(_CANONICAL_CHR_ORDER)}


def _validate_reference_sort_order(chr_col: pd.Series, bp_col: np.ndarray, path: Path) -> None:
    chr_rank = chr_col.map(_CHR_RANK)
    bad_chr = chr_rank.isna()
    if bad_chr.any():
        idx = int(bad_chr.idxmax())
        raise ValueError(
            f"{path}:{idx + 1}: chr must use canonical labels 1-22 or X"
        )

    chr_rank_vals = chr_rank.to_numpy(dtype=np.int64)
    bp_vals = np.asarray(bp_col, dtype=np.int64)
    bad_order = (chr_rank_vals[1:] < chr_rank_vals[:-1]) | (
        (chr_rank_vals[1:] == chr_rank_vals[:-1]) & (bp_vals[1:] < bp_vals[:-1])
    )
    if bad_order.any():
        idx = int(np.flatnonzero(bad_order)[0]) + 1
        raise ValueError(
            f"{path}:{idx + 1}: rows must be sorted by canonical chromosome (1-22, X) and ascending bp"
        )


def _parse_bim(path: Path) -> list:
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

    bad_allele = a1_col.eq("") | a2_col.eq("")
    if bad_allele.any():
        idx = int(bad_allele.idxmax())
        raise ValueError(f"{path}:{idx + 1}: a1 and a2 must be non-empty")

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
    _validate_reference_sort_order(chr_col, bp_int.to_numpy(dtype=np.int64), path)

    # cm is validated for BIM schema compatibility but intentionally ignored
    # in-memory per spec.
    return list(
        zip(
            chr_col.tolist(),
            snp_col.tolist(),
            bp_int.tolist(),
            a1_col.tolist(),
            a2_col.tolist(),
        )
    )


def _checksum(rows) -> str:
    text = "".join(f"{r[0]}:{r[2]}:{r[3]}:{r[4]}\n" for r in rows)
    return hashlib.md5(text.encode()).hexdigest()


def _chr_sort_key(label: str):
    if label not in _CHR_RANK:
        raise ValueError(f"Unsupported shard label: {label!r}; expected canonical labels 1-22 or X")
    return _CHR_RANK[label]


class ReferenceShard:
    def __init__(self, label: str, rows, *, _checksum_override: str | None = None):
        self._label = label
        self._chr = np.array([r[0] for r in rows], dtype=object)
        self._snp = np.array([r[1] for r in rows], dtype=object)
        self._bp = np.array([r[2] for r in rows], dtype=np.int64)
        self._a1 = np.array([r[3] for r in rows], dtype=object)
        self._a2 = np.array([r[4] for r in rows], dtype=object)
        self._checksum = (
            _checksum_override if _checksum_override is not None else _checksum(rows)
        )

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


def load_reference(path) -> ReferencePanel:
    path_str = str(path)
    path_obj = Path(path_str)

    if "@" in path_str:
        parent = path_obj.parent
        name = path_obj.name
        at_idx = name.index("@")
        name_prefix = name[:at_idx]
        name_suffix = name[at_idx + 1:]

        candidates = []
        for f in sorted(parent.iterdir()):
            if not f.is_file():
                continue
            n = f.name
            if not n.startswith(name_prefix):
                continue
            if name_suffix and not n.endswith(name_suffix):
                continue
            inner = n[len(name_prefix): -len(name_suffix) if name_suffix else None]
            if inner:
                candidates.append((inner, f))

        try:
            candidates.sort(key=lambda x: _chr_sort_key(x[0]))
        except ValueError as exc:
            raise ValueError(f"{path_str}: {exc}") from None

        if not candidates:
            raise FileNotFoundError(
                f"No BIM shards found matching template: {path_str}"
            )

        shards = [
            ReferenceShard(label, _parse_bim(bim_path))
            for label, bim_path in candidates
        ]
        return ReferencePanel(shards)

    rows = _parse_bim(path_obj)

    by_chr: dict = {}
    for r in rows:
        c = r[0]
        if c not in by_chr:
            by_chr[c] = []
        by_chr[c].append(r)

    return ReferencePanel([ReferenceShard(c, by_chr[c]) for c in by_chr])


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
        arrays[p + "chr"] = np.asarray([str(v) for v in s.chr], dtype=str)
        arrays[p + "snp"] = np.asarray([str(v) for v in s.snp], dtype=str)
        arrays[p + "bp"]  = s.bp
        arrays[p + "a1"]  = np.asarray([str(v) for v in s.a1], dtype=str)
        arrays[p + "a2"]  = np.asarray([str(v) for v in s.a2], dtype=str)
    np.savez_compressed(path, **arrays)


def load_reference_cache(path) -> ReferencePanel:
    data = np.load(path, allow_pickle=False)
    meta = json.loads(bytes(data["_meta"]).decode())
    schema = meta.get("schema")
    if schema != _CACHE_SCHEMA:
        raise ValueError(f"Unsupported reference cache schema: {schema!r}")
    shards = []
    for i, (label, checksum) in enumerate(
        zip(meta["shard_labels"], meta["shard_checksums"])
    ):
        p = f"s{i}_"
        shards.append(
            ReferenceShard._from_arrays(
                label,
                data[p + "chr"],
                data[p + "snp"],
                data[p + "bp"],
                data[p + "a1"],
                data[p + "a2"],
                checksum,
            )
        )
    return ReferencePanel(shards)
