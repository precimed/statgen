import json
import io
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError, ParserError
from scipy import sparse

from ._utils import validate_requested_shards

_CACHE_SCHEMA = "annotations_cache/0.1"


def _is_ignored_bed_line(line: str) -> bool:
    return (
        line == ""
        or line.startswith("#")
        or line.startswith("track ")
        or line.startswith("browser ")
    )


def _read_bed_payload(path: Path) -> str:
    kept = []
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.rstrip("\r\n")
            if _is_ignored_bed_line(line):
                continue
            kept.append(line)
    if not kept:
        raise ValueError(f"{path}: BED file is empty")
    return "\n".join(kept) + "\n"


def _as_sparse_binary(mat) -> sparse.csr_matrix:
    if sparse.issparse(mat):
        csr = mat.tocsr(copy=True)
    else:
        arr = np.asarray(mat)
        if arr.ndim != 2:
            raise ValueError("annomat must be a 2D matrix")
        csr = sparse.csr_matrix(arr)

    # Guard against CSR matrices carrying explicit stored zeros in .data.
    # These represent logical zeros and must not be flipped to ones.
    csr.eliminate_zeros()

    if csr.data.size:
        data = np.asarray(csr.data)
        bad = data != 1
        if bad.any():
            idx = int(np.flatnonzero(bad)[0])
            raise ValueError(f"annomat contains non-binary value: {data[idx]!r}")

    csr = csr.astype(np.uint8)
    csr.data[:] = 1
    return csr


def _coerce_annovec(annovec, n: int) -> np.ndarray:
    arr = np.asarray(annovec)
    if arr.ndim > 2 or (arr.ndim == 2 and 1 not in arr.shape):
        raise ValueError(f"annovec must be a vector with length {n}")
    arr = arr.reshape(-1)
    if arr.size != n:
        raise ValueError(f"annovec length mismatch: expected {n}, got {arr.size}")
    bad = (arr != 0) & (arr != 1)
    if bad.any():
        idx = int(np.flatnonzero(bad)[0])
        raise ValueError(f"annovec[{idx}] must be binary (0/1)")
    return arr.astype(np.uint8)


def _coerce_annonames(annonames) -> list[str]:
    if isinstance(annonames, str):
        raise ValueError("annonames must be a non-empty list of unique strings")
    names = [str(x) for x in list(annonames)]
    if not names:
        raise ValueError("annonames must be a non-empty list of unique strings")
    if any(n == "" for n in names):
        raise ValueError("annonames must not contain empty strings")
    if len(set(names)) != len(names):
        raise ValueError("annonames must be unique")
    return names


def _annotation_name_from_path(path: Path) -> str:
    name = path.stem
    if not name:
        raise ValueError(f"invalid annotation filename: {path}")
    return name


def _parse_bed(path: Path) -> dict[str, np.ndarray]:
    payload = _read_bed_payload(path)
    try:
        df = pd.read_csv(
            io.StringIO(payload),
            sep=r"\s+",
            header=None,
            dtype=str,
            engine="python",
            keep_default_na=False,
            na_filter=False,
        )
    except EmptyDataError as exc:
        raise ValueError(f"{path}: BED file is empty") from exc
    except ParserError as exc:
        raise ValueError(f"{path}: malformed BED") from exc

    if df.shape[0] == 0:
        raise ValueError(f"{path}: BED file is empty")
    if df.shape[1] < 3:
        raise ValueError(f"{path}: BED must have at least 3 whitespace-separated columns")

    chr_col = df.iloc[:, 0]
    start_raw = df.iloc[:, 1]
    end_raw = df.iloc[:, 2]

    bad_chr = chr_col.eq("")
    if bad_chr.any():
        idx = int(bad_chr.idxmax())
        raise ValueError(f"{path}: row {idx + 1}: chromosome label must be non-empty")

    start_num = pd.to_numeric(start_raw, errors="coerce")
    end_num = pd.to_numeric(end_raw, errors="coerce")

    bad_start = start_num.isna() | (np.floor(start_num) != start_num) | (start_num < 0)
    if bad_start.any():
        idx = int(bad_start.idxmax())
        raise ValueError(f"{path}: row {idx + 1}: BED start must be a non-negative integer")

    bad_end = end_num.isna() | (np.floor(end_num) != end_num) | (end_num < 0)
    if bad_end.any():
        idx = int(bad_end.idxmax())
        raise ValueError(f"{path}: row {idx + 1}: BED end must be a non-negative integer")

    start_int = start_num.astype(np.int64)
    end_int = end_num.astype(np.int64)

    bad_len = end_int < start_int
    if bad_len.any():
        idx = int(bad_len.idxmax())
        raise ValueError(f"{path}: row {idx + 1}: BED interval end must be >= start")

    intervals_by_chr: dict[str, np.ndarray] = {}
    for chr_label, grp in pd.DataFrame(
        {"chr": chr_col.to_numpy(dtype=object), "start": start_int, "end": end_int}
    ).groupby("chr", sort=False):
        starts = grp["start"].to_numpy(dtype=np.int64)
        ends = grp["end"].to_numpy(dtype=np.int64)
        intervals_by_chr[str(chr_label)] = _merge_intervals(starts, ends)

    return intervals_by_chr


def _merge_intervals(starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    if starts.size == 0:
        return np.empty((0, 2), dtype=np.int64)

    order = np.lexsort((ends, starts))
    s = starts[order]
    e = ends[order]

    merged_start = [int(s[0])]
    merged_end = [int(e[0])]
    for i in range(1, s.size):
        cur_s = int(s[i])
        cur_e = int(e[i])
        if cur_s <= merged_end[-1]:
            if cur_e > merged_end[-1]:
                merged_end[-1] = cur_e
        else:
            merged_start.append(cur_s)
            merged_end.append(cur_e)

    return np.column_stack(
        [np.asarray(merged_start, dtype=np.int64), np.asarray(merged_end, dtype=np.int64)]
    )


def _paint_mask(bp: np.ndarray, intervals: np.ndarray) -> np.ndarray:
    n = bp.size
    if intervals.size == 0:
        return np.zeros(n, dtype=np.uint8)

    pos0 = np.asarray(bp, dtype=np.int64) - 1
    starts = intervals[:, 0]
    ends = intervals[:, 1]

    idx = np.searchsorted(starts, pos0, side="right") - 1
    valid = idx >= 0
    out = np.zeros(n, dtype=np.uint8)
    if valid.any():
        valid_idx = idx[valid]
        inside = pos0[valid] < ends[valid_idx]
        out[np.flatnonzero(valid)[inside]] = 1
    return out


def _paint_annotations(bed_paths: list[Path], reference) -> tuple[sparse.csr_matrix, list[str]]:
    annonames = [_annotation_name_from_path(p) for p in bed_paths]
    if len(set(annonames)) != len(annonames):
        raise ValueError("duplicate annotation names derived from BED basenames")

    n = int(reference.num_snp)
    k = len(bed_paths)
    if n == 0:
        return sparse.csr_matrix((0, k), dtype=np.uint8), annonames

    columns = []
    for bed_path in bed_paths:
        intervals_by_chr = _parse_bed(bed_path)
        mask = np.zeros(n, dtype=np.uint8)
        for ref_shard, off in zip(reference.shards, reference.shard_offsets):
            start = int(off["start0"])
            stop = int(off["stop0"])
            intervals = intervals_by_chr.get(ref_shard.label)
            if intervals is None:
                continue
            mask[start:stop] = _paint_mask(ref_shard.bp, intervals)
        columns.append(mask)

    dense = np.column_stack(columns).astype(np.uint8)
    return sparse.csr_matrix(dense), annonames


class AnnotationShard:
    def __init__(self, label: str, checksum: str, annomat):
        self._label = str(label)
        self._checksum = str(checksum)
        self._annomat = _as_sparse_binary(annomat)

    @property
    def label(self) -> str:
        return self._label

    @property
    def checksum(self) -> str:
        return self._checksum

    @property
    def num_snp(self) -> int:
        return int(self._annomat.shape[0])

    @property
    def num_annot(self) -> int:
        return int(self._annomat.shape[1])

    @property
    def annomat(self) -> sparse.csr_matrix:
        return self._annomat

    @classmethod
    def _from_arrays(cls, label, checksum, annomat):
        return cls(label=label, checksum=checksum, annomat=annomat)


class AnnotationPanel:
    def __init__(self, shards: list[AnnotationShard], annonames):
        self._shards = list(shards)
        self._annonames = np.asarray(_coerce_annonames(annonames), dtype=object)
        self._num_snp = int(sum(s.num_snp for s in self._shards))

        self._shard_offsets = []
        pos = 0
        for s in self._shards:
            if s.num_annot != self._annonames.size:
                raise ValueError("all shards must share the same annotation columns")
            self._shard_offsets.append(
                {"shard_label": s.label, "start0": pos, "stop0": pos + s.num_snp}
            )
            pos += s.num_snp

    @property
    def shards(self) -> list[AnnotationShard]:
        return list(self._shards)

    @property
    def annonames(self) -> np.ndarray:
        return self._annonames.copy()

    @property
    def num_snp(self) -> int:
        return self._num_snp

    @property
    def num_annot(self) -> int:
        return int(self._annonames.size)

    @property
    def shard_offsets(self) -> list:
        return list(self._shard_offsets)

    @property
    def annomat(self) -> sparse.csr_matrix:
        if not self._shards:
            return sparse.csr_matrix((0, self.num_annot), dtype=np.uint8)
        return sparse.vstack([s.annomat for s in self._shards], format="csr")

    def select_shards(self, shards) -> "AnnotationPanel":
        available = [s.label for s in self._shards]
        selected = validate_requested_shards(shards, available, "AnnotationPanel.select_shards")
        by_label = {s.label: s for s in self._shards}
        return AnnotationPanel([by_label[label] for label in selected], self._annonames)

    def select_annotations(self, names) -> "AnnotationPanel":
        if isinstance(names, str):
            raise ValueError("names must be a non-empty list of unique annotation names")
        names_list = [str(x) for x in list(names)]
        if not names_list:
            raise ValueError("names must be a non-empty list of unique annotation names")
        if len(set(names_list)) != len(names_list):
            raise ValueError("names must be unique")

        idx_map = {name: i for i, name in enumerate(self._annonames.tolist())}
        missing = [name for name in names_list if name not in idx_map]
        if missing:
            raise ValueError(f"unknown annotation name(s): {', '.join(missing)}")

        idx = np.asarray([idx_map[name] for name in names_list], dtype=np.int64)
        out_shards = [
            AnnotationShard._from_arrays(s.label, s.checksum, s.annomat[:, idx])
            for s in self._shards
        ]
        return AnnotationPanel(out_shards, names_list)

    def union_annotations(self, other, mode: str = "by_name") -> "AnnotationPanel":
        if mode != "by_name":
            raise ValueError("union_annotations supports only mode='by_name'")

        lhs_names = self._annonames.tolist()
        rhs_names = np.asarray(getattr(other, "annonames", []), dtype=object).tolist()
        overlap = sorted(set(lhs_names).intersection(rhs_names))
        if overlap:
            raise ValueError(f"annotation name collision(s): {', '.join(overlap)}")

        other_shards = list(getattr(other, "shards", []))
        if len(other_shards) != len(self._shards):
            raise ValueError("union_annotations requires matching shard structure")

        out_shards = []
        for a, b in zip(self._shards, other_shards):
            if a.label != getattr(b, "label", None):
                raise ValueError("union_annotations requires matching shard labels")
            if a.num_snp != getattr(b, "num_snp", None):
                raise ValueError("union_annotations requires matching shard row counts")
            b_checksum = getattr(b, "checksum", None)
            if b_checksum is not None and b_checksum != a.checksum:
                raise ValueError("union_annotations requires checksum-compatible reference alignment")

            b_anno = getattr(b, "annomat", None)
            if b_anno is None:
                raise ValueError("union_annotations requires other shards to expose annomat")
            union_mat = sparse.hstack([a.annomat, _as_sparse_binary(b_anno)], format="csr")
            out_shards.append(AnnotationShard._from_arrays(a.label, a.checksum, union_mat))

        return AnnotationPanel(out_shards, lhs_names + rhs_names)


def create_annotations(reference, annomat, annonames) -> AnnotationPanel:
    names = _coerce_annonames(annonames)
    n = int(reference.num_snp)
    mat = _as_sparse_binary(annomat)
    if mat.shape != (n, len(names)):
        raise ValueError(
            f"annomat shape mismatch: expected ({n}, {len(names)}), got {mat.shape}"
        )

    out_shards = []
    for ref_shard, off in zip(reference.shards, reference.shard_offsets):
        start = int(off["start0"])
        stop = int(off["stop0"])
        out_shards.append(
            AnnotationShard(
                label=ref_shard.label,
                checksum=ref_shard.checksum,
                annomat=mat[start:stop, :],
            )
        )
    return AnnotationPanel(out_shards, names)


def create_annotation(reference, annovec, annoname) -> AnnotationPanel:
    name = str(annoname)
    if not name:
        raise ValueError("annoname must be non-empty")
    vec = _coerce_annovec(annovec, int(reference.num_snp))
    mat = sparse.csr_matrix(vec.reshape(-1, 1))
    return create_annotations(reference, mat, [name])


def load_annotations(bed_paths, reference) -> AnnotationPanel:
    if isinstance(bed_paths, (str, Path)):
        paths = [Path(bed_paths)]
    else:
        paths = [Path(p) for p in list(bed_paths)]
        if not paths:
            raise ValueError("bed_paths must be a non-empty list of BED files")
    for p in paths:
        if not p.is_file():
            raise FileNotFoundError(f"BED file not found: {p}")

    annomat, annonames = _paint_annotations(paths, reference)
    return create_annotations(reference, annomat, annonames)


def save_annotations_cache(panel: AnnotationPanel, path) -> None:
    meta = {
        "schema": _CACHE_SCHEMA,
        "shard_labels": [s.label for s in panel.shards],
        "shard_checksums": [s.checksum for s in panel.shards],
        "annonames": panel.annonames.tolist(),
    }
    arrays = {
        "_meta": np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8),
    }

    for i, shard in enumerate(panel.shards):
        p = f"s{i}_"
        mat = shard.annomat.tocsr()
        arrays[p + "data"] = mat.data.astype(np.uint8)
        arrays[p + "indices"] = mat.indices.astype(np.int32)
        arrays[p + "indptr"] = mat.indptr.astype(np.int32)
        arrays[p + "shape"] = np.asarray(mat.shape, dtype=np.int64)

    np.savez_compressed(path, **arrays)


def load_annotations_cache(path, shards=None) -> AnnotationPanel:
    with np.load(path, allow_pickle=False) as data:
        meta = json.loads(bytes(data["_meta"]).decode())
        schema = meta.get("schema")
        if schema != _CACHE_SCHEMA:
            raise ValueError(f"Unsupported annotations cache schema: {schema!r}")

        labels = list(meta.get("shard_labels", []))
        checksums = list(meta.get("shard_checksums", []))
        if len(labels) != len(checksums):
            raise ValueError(
                "Invalid annotations cache: shard_labels and shard_checksums length mismatch"
            )

        annonames = _coerce_annonames(meta.get("annonames", []))
        selected = validate_requested_shards(shards, labels, "load_annotations_cache")
        label_to_index = {label: i for i, label in enumerate(labels)}

        shard_objs = []
        for label in selected:
            i = label_to_index[label]
            p = f"s{i}_"
            mat = sparse.csr_matrix(
                (
                    np.asarray(data[p + "data"], dtype=np.uint8),
                    np.asarray(data[p + "indices"], dtype=np.int32),
                    np.asarray(data[p + "indptr"], dtype=np.int32),
                ),
                shape=tuple(np.asarray(data[p + "shape"], dtype=np.int64).tolist()),
            )
            shard_objs.append(
                AnnotationShard._from_arrays(
                    label=label,
                    checksum=checksums[i],
                    annomat=mat,
                )
            )

    return AnnotationPanel(shard_objs, annonames)
