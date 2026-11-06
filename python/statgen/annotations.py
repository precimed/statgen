from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError, ParserError
from scipy import sparse

from ._annotation_painting import (
    IGNORED_CHR,
    binary_intervals_by_chr_from_arrays,
    paint_binary_column,
    paint_mask,
    paint_numeric_sparse,
    validate_numeric_non_overlapping,
)
from .reference import _raise_if_reference_shards_incompatible
from ._utils import CANONICAL_CHR_ORDER, validate_requested_shards

_CACHE_SCHEMA = "annotations_cache/0.2"
_OLD_CACHE_SCHEMA = "annotations_cache/0.1"
_ALLOWED_CHR = set(CANONICAL_CHR_ORDER) | IGNORED_CHR


def _as_sparse_numeric(mat, name: str = "annomat") -> sparse.csr_matrix:
    if sparse.issparse(mat):
        csr = mat.tocsr(copy=True)
    else:
        arr = np.asarray(mat)
        if arr.ndim != 2:
            raise ValueError(f"{name} must be a 2D matrix")
        csr = sparse.csr_matrix(arr)

    csr = csr.astype(np.float64)
    csr.eliminate_zeros()
    if csr.data.size and not np.isfinite(csr.data).all():
        idx = int(np.flatnonzero(~np.isfinite(csr.data))[0])
        raise ValueError(f"{name} contains non-finite value: {csr.data[idx]!r}")
    return csr


def _coerce_annovec(annovec, n: int) -> np.ndarray:
    arr = np.asarray(annovec, dtype=np.float64)
    if arr.ndim > 2 or (arr.ndim == 2 and 1 not in arr.shape):
        raise ValueError(f"annovec must be a vector with length {n}")
    arr = arr.reshape(-1)
    if arr.size != n:
        raise ValueError(f"annovec length mismatch: expected {n}, got {arr.size}")
    if not np.isfinite(arr).all():
        idx = int(np.flatnonzero(~np.isfinite(arr))[0])
        raise ValueError(f"annovec[{idx}] must be finite numeric")
    return arr


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


def _coerce_string_vector(value, expected_len: int, name: str) -> np.ndarray:
    if isinstance(value, str):
        raise ValueError(f"{name} must be a string vector with length {expected_len}")
    values = [str(x) for x in list(value)]
    if len(values) != expected_len:
        raise ValueError(f"{name} length mismatch: expected {expected_len}, got {len(values)}")
    return np.asarray(values, dtype=object)


def _coerce_is_binary(is_binary, mat: sparse.csr_matrix) -> np.ndarray:
    n_annot = int(mat.shape[1])
    if is_binary is None:
        out = np.zeros(n_annot, dtype=bool)
        csc = mat.tocsc(copy=False)
        for j in range(n_annot):
            start = int(csc.indptr[j])
            stop = int(csc.indptr[j + 1])
            data = csc.data[start:stop]
            out[j] = bool(data.size == 0 or np.all(data == 1.0))
        return out

    arr = np.asarray(is_binary)
    if arr.ndim == 0:
        raise ValueError(f"is_binary must be a logical vector with length {n_annot}")
    arr = arr.reshape(-1)
    if arr.size != n_annot:
        raise ValueError(f"is_binary length mismatch: expected {n_annot}, got {arr.size}")
    if arr.dtype != np.bool_:
        bad = ~np.isin(arr, [False, True, 0, 1])
        if bad.any():
            idx = int(np.flatnonzero(bad)[0])
            raise ValueError(f"is_binary[{idx}] must be logical")
    return arr.astype(bool)


def _coerce_is_binary_vector(is_binary, n_annot: int) -> np.ndarray:
    if is_binary is None:
        raise ValueError("is_binary must be supplied")
    arr = np.asarray(is_binary)
    if arr.ndim == 0:
        raise ValueError(f"is_binary must be a logical vector with length {n_annot}")
    arr = arr.reshape(-1)
    if arr.size != n_annot:
        raise ValueError(f"is_binary length mismatch: expected {n_annot}, got {arr.size}")
    if arr.dtype != np.bool_:
        bad = ~np.isin(arr, [False, True, 0, 1])
        if bad.any():
            idx = int(np.flatnonzero(bad)[0])
            raise ValueError(f"is_binary[{idx}] must be logical")
    return arr.astype(bool)


def _validate_declared_binary(mat: sparse.csr_matrix, is_binary: np.ndarray, name: str) -> None:
    if not np.any(is_binary):
        return
    data = mat[:, is_binary].data
    if data.size:
        bad = data != 1.0
        if bad.any():
            idx = int(np.flatnonzero(bad)[0])
            raise ValueError(
                f"{name} declared binary must be binary (0/1); found non-binary value {data[idx]!r}"
            )


def _annotation_name_from_path(path: Path) -> str:
    name = path.stem
    if not name:
        raise ValueError(f"invalid annotation filename: {path}")
    return name


def _generated_metadata(
    path: Path,
    has_header: bool,
    num_source_intervals: int,
    *,
    value_column=None,
    group_column=None,
    group_value=None,
) -> str:
    payload = {
        "source_file": str(path),
        "source_file_has_header": 1 if has_header else 0,
        "num_source_intervals": int(num_source_intervals),
    }
    if value_column is not None:
        payload["value_column"] = value_column
    if group_column is not None:
        payload["group_column"] = group_column
    if group_value is not None:
        payload["group_value"] = group_value
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


def _read_sidecar_exact(path) -> str:
    return Path(path).read_text(encoding="utf-8")


def _read_column_metadata_sidecar(path, num_columns: int) -> list[str]:
    text = _read_sidecar_exact(path)
    lines = text.splitlines()
    if len(lines) != num_columns:
        raise ValueError(
            f"{path}: annotation metadata sidecar line count mismatch: "
            f"expected {num_columns}, got {len(lines)}"
        )
    empty = [i for i, line in enumerate(lines) if line == ""]
    if empty:
        raise ValueError(f"{path}: annotation metadata sidecar line {empty[0] + 1} is empty")
    return lines


def _read_annotation_table(path: Path, has_header: bool) -> tuple[pd.DataFrame, list[str] | None, int]:
    if not path.is_file():
        raise FileNotFoundError(f"annotation file not found: {path}")

    try:
        df = pd.read_csv(
            path,
            sep="\t",
            header=0 if has_header else None,
            dtype=str,
            keep_default_na=False,
            na_filter=False,
            comment="#",
            skip_blank_lines=True,
        )
    except EmptyDataError as exc:
        raise ValueError(f"{path}: BED file is empty") from exc
    except ParserError as exc:
        raise ValueError(f"{path}: malformed annotation input") from exc

    if df.shape[0] == 0:
        raise ValueError(f"{path}: BED file is empty")

    if df.shape[1] < 3:
        raise ValueError(f"{path}: BED must have at least 3 tab-separated columns")
    header_fields = [str(x) for x in df.columns] if has_header else None
    if header_fields is not None:
        if any(x == "" for x in header_fields):
            raise ValueError(f"{path}: header names must be non-empty")
        if len(set(header_fields)) != len(header_fields):
            raise ValueError(f"{path}: header names must be unique")

    return df, header_fields, 1 if has_header else 0


def _validate_interval_columns(
    df: pd.DataFrame,
    path: Path,
    row_base0: int = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    chr_col = df.iloc[:, 0].astype(str)
    bad_chr = chr_col.eq("")
    if bad_chr.any():
        idx = int(bad_chr.idxmax())
        raise ValueError(f"{path}: row {row_base0 + idx + 1}: chromosome label must be non-empty")
    chr_values = chr_col.to_numpy(dtype=object)

    chr_series = pd.Series(chr_values, dtype="string")
    chr_style = chr_series.str.startswith("chr", na=False)
    if chr_style.any():
        idx = int(chr_style.idxmax())
        raise ValueError(f"{path}: row {row_base0 + idx + 1}: chr-style labels (e.g., chr1/chrX) are not allowed")
    unsupported = ~chr_series.isin(_ALLOWED_CHR)
    if unsupported.any():
        idx = int(unsupported.idxmax())
        raise ValueError(
            f"{path}: row {row_base0 + idx + 1}: unsupported chr label {chr_values[idx]!r}; expected 1-22, X (Y/MT are ignored)"
        )

    start_raw = df.iloc[:, 1]
    end_raw = df.iloc[:, 2]
    start_num = pd.to_numeric(start_raw, errors="coerce")
    end_num = pd.to_numeric(end_raw, errors="coerce")

    bad_start = start_num.isna() | (np.floor(start_num) != start_num) | (start_num < 0)
    if bad_start.any():
        idx = int(bad_start.idxmax())
        raise ValueError(f"{path}: row {row_base0 + idx + 1}: BED start must be a non-negative integer")

    bad_end = end_num.isna() | (np.floor(end_num) != end_num) | (end_num < 0)
    if bad_end.any():
        idx = int(bad_end.idxmax())
        raise ValueError(f"{path}: row {row_base0 + idx + 1}: BED end must be a non-negative integer")

    start_int = start_num.astype(np.int64).to_numpy()
    end_int = end_num.astype(np.int64).to_numpy()
    bad_len = end_int < start_int
    if bad_len.any():
        idx = int(np.flatnonzero(bad_len)[0])
        raise ValueError(f"{path}: row {row_base0 + idx + 1}: BED interval end must be >= start")

    return chr_values, start_int, end_int


def _binary_intervals_by_chr(path: Path) -> tuple[dict[str, np.ndarray], int]:
    df, _, row_base0 = _read_annotation_table(path, has_header=False)
    if df.shape[1] < 3:
        raise ValueError(f"{path}: BED must have at least 3 tab-separated columns")
    chr_values, starts, ends = _validate_interval_columns(df, path, row_base0=row_base0)
    return binary_intervals_by_chr_from_arrays(chr_values, starts, ends), int(df.shape[0])


def _paint_binary_annotations(
    bed_paths: list[Path], reference
) -> tuple[sparse.csr_matrix, list[str], list[int]]:
    annonames = [_annotation_name_from_path(p) for p in bed_paths]
    if len(set(annonames)) != len(annonames):
        raise ValueError("duplicate annotation names derived from BED basenames")

    n = int(reference.num_snp)
    k = len(bed_paths)
    if n == 0:
        counts = []
        for bed_path in bed_paths:
            _, count = _binary_intervals_by_chr(bed_path)
            counts.append(count)
        return sparse.csr_matrix((0, k), dtype=np.float64), annonames, counts

    columns = []
    counts = []
    for bed_path in bed_paths:
        intervals_by_chr, count = _binary_intervals_by_chr(bed_path)
        columns.append(paint_binary_column(intervals_by_chr, reference))
        counts.append(count)

    return sparse.hstack(columns, format="csr"), annonames, counts


def _selector_list(selectors, name: str) -> list | None:
    if selectors is None:
        return None
    if isinstance(selectors, (str, int, np.integer)):
        values = [selectors]
    else:
        values = list(selectors)
    if not values:
        raise ValueError(f"{name} must be a non-empty list")
    return values


def _normalize_column_selectors(selectors, header_fields, num_columns: int, path: Path, name: str) -> tuple[list[int], list] | None:
    selector_values = _selector_list(selectors, name)
    if selector_values is None:
        return None

    out = []
    metadata_values = []
    singular = name[:-1] if name.endswith("s") else name
    for selector in selector_values:
        if isinstance(selector, str):
            if header_fields is None:
                raise ValueError(f"named {name} are invalid when has_header=False")
            if selector not in header_fields:
                raise ValueError(f"{path}: unknown {singular} name {selector!r}")
            col0 = header_fields.index(selector)
            metadata_value = selector
        elif isinstance(selector, (int, np.integer)) and not isinstance(selector, bool):
            col0 = int(selector)
            metadata_value = int(selector)
        else:
            raise ValueError(f"{name} must contain column names or integer indices")
        if col0 < 3 or col0 >= num_columns:
            raise ValueError(
                f"{path}: {singular} {selector!r} is out of range; selected columns must be physical columns 4 or later"
            )
        out.append(col0)
        metadata_values.append(metadata_value)

    if len(set(out)) != len(out):
        raise ValueError(f"{name} must not contain duplicates")
    return out, metadata_values


def _normalize_value_columns(value_columns, header_fields, num_columns: int, path: Path) -> tuple[list[int], list] | None:
    return _normalize_column_selectors(value_columns, header_fields, num_columns, path, "value_columns")


def _normalize_group_column(group_column, header_fields, num_columns: int, path: Path) -> tuple[int, object] | None:
    normalized = _normalize_column_selectors(group_column, header_fields, num_columns, path, "group_column")
    if normalized is None:
        return None
    columns, metadata_values = normalized
    if len(columns) != 1:
        raise ValueError("group_column must identify exactly one source column")
    return columns[0], metadata_values[0]


def _default_annotation_names(
    path: Path,
    num_columns: int,
    header_fields: list[str] | None,
    value_columns0: list[int] | None,
) -> list[str]:
    if value_columns0 is None:
        if num_columns == 3:
            return [_annotation_name_from_path(path)]
        if num_columns == 4:
            if header_fields is None:
                return [_annotation_name_from_path(path)]
            return [str(header_fields[3])]
        raise ValueError(f"{path}: input with five or more columns requires explicit value_columns or group_column")

    if header_fields is None and num_columns >= 5:
        raise ValueError(f"{path}: headerless input with five or more columns requires annotation_names")
    if header_fields is not None:
        return [str(header_fields[i]) for i in value_columns0]
    return [_annotation_name_from_path(path)]


def _numeric_values(
    df: pd.DataFrame,
    value_columns0: list[int],
    path: Path,
    row_base0: int = 0,
) -> np.ndarray:
    selected = df.iloc[:, value_columns0]
    numeric = selected.apply(pd.to_numeric, errors="coerce")
    values = numeric.to_numpy(dtype=np.float64)
    finite = np.isfinite(values)
    if not finite.all():
        row0, col_pos = np.argwhere(~finite)[0]
        col0 = value_columns0[int(col_pos)]
        raise ValueError(
            f"{path}: row {row_base0 + int(row0) + 1}: annotation value column {col0} must be finite numeric"
        )
    return values


class AnnotationShard:
    def __init__(self, label: str, reference_checksum: str, annomat):
        self._label = str(label)
        self._reference_checksum = str(reference_checksum)
        self._annomat = _as_sparse_numeric(annomat)

    @property
    def label(self) -> str:
        return self._label

    @property
    def reference_checksum(self) -> str:
        return self._reference_checksum

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
    def _from_arrays(cls, label, reference_checksum, annomat):
        return cls(label=label, reference_checksum=reference_checksum, annomat=annomat)


class AnnotationPanel:
    def __init__(
        self,
        shards: list[AnnotationShard],
        annonames,
        is_binary=None,
        annotation_metadata=None,
    ):
        self._shards = list(shards)
        if not self._shards:
            raise ValueError("AnnotationPanel requires at least one shard")
        self._annonames = np.asarray(_coerce_annonames(annonames), dtype=object)
        self._num_snp = int(sum(s.num_snp for s in self._shards))

        for s in self._shards:
            if s.num_annot != self._annonames.size:
                raise ValueError("all shards must share the same annotation columns")

        self._is_binary = _coerce_is_binary_vector(is_binary, self._annonames.size)

        if annotation_metadata is None:
            self._annotation_metadata = np.asarray([""] * self._annonames.size, dtype=object)
        else:
            self._annotation_metadata = _coerce_string_vector(
                annotation_metadata, self._annonames.size, "annotation_metadata"
            )

        self._shard_offsets = []
        pos = 0
        for s in self._shards:
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
    def is_binary(self) -> np.ndarray:
        return self._is_binary.copy()

    @property
    def annotation_metadata(self) -> np.ndarray:
        return self._annotation_metadata.copy()

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
            return sparse.csr_matrix((0, self.num_annot), dtype=np.float64)
        return sparse.vstack([s.annomat for s in self._shards], format="csr")

    def select_shards(self, shards) -> "AnnotationPanel":
        available = [s.label for s in self._shards]
        selected = validate_requested_shards(shards, available, "AnnotationPanel.select_shards")
        by_label = {s.label: s for s in self._shards}
        return AnnotationPanel(
            [by_label[label] for label in selected],
            self._annonames,
            is_binary=self._is_binary,
            annotation_metadata=self._annotation_metadata,
        )

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
            AnnotationShard._from_arrays(s.label, s.reference_checksum, s.annomat[:, idx])
            for s in self._shards
        ]
        return AnnotationPanel(
            out_shards,
            names_list,
            is_binary=self._is_binary[idx],
            annotation_metadata=self._annotation_metadata[idx],
        )

    def union_annotations(self, other, mode: str = "by_name") -> "AnnotationPanel":
        if mode != "by_name":
            raise ValueError("union_annotations supports only mode='by_name'")
        if not isinstance(other, AnnotationPanel):
            raise ValueError("union_annotations requires another AnnotationPanel")

        lhs_names = self._annonames.tolist()
        rhs_names = other._annonames.tolist()
        overlap = sorted(set(lhs_names).intersection(rhs_names))
        if overlap:
            raise ValueError(f"annotation name collision(s): {', '.join(overlap)}")

        _raise_if_reference_shards_incompatible(
            self._shards,
            other._shards,
            where="union_annotations",
            require_checksum=True,
        )

        out_shards = []
        for a, b in zip(self._shards, other._shards):
            union_mat = sparse.hstack([a.annomat, b.annomat], format="csr")
            out_shards.append(AnnotationShard._from_arrays(a.label, a.reference_checksum, union_mat))

        return AnnotationPanel(
            out_shards,
            lhs_names + rhs_names,
            is_binary=np.concatenate([self._is_binary, other._is_binary]),
            annotation_metadata=np.concatenate([self._annotation_metadata, other._annotation_metadata]),
        )

    def save_cache(self, path) -> None:
        save_annotations_cache(self, path)


def create_annotations(
    reference,
    annotation_matrix,
    annotation_names,
    is_binary=None,
    annotation_metadata=None,
) -> AnnotationPanel:
    names = _coerce_annonames(annotation_names)
    n = int(reference.num_snp)
    mat = _as_sparse_numeric(annotation_matrix, "annotation_matrix")
    if mat.shape != (n, len(names)):
        raise ValueError(
            f"annotation_matrix shape mismatch: expected ({n}, {len(names)}), got {mat.shape}"
        )

    binary = _coerce_is_binary(is_binary, mat)
    _validate_declared_binary(mat, binary, "annotation_matrix")
    metadata = (
        np.asarray([""] * len(names), dtype=object)
        if annotation_metadata is None
        else _coerce_string_vector(annotation_metadata, len(names), "annotation_metadata")
    )

    out_shards = []
    for ref_shard, off in zip(reference.shards, reference.shard_offsets):
        start = int(off["start0"])
        stop = int(off["stop0"])
        out_shards.append(
            AnnotationShard(
                label=ref_shard.label,
                reference_checksum=ref_shard.checksum,
                annomat=mat[start:stop, :],
            )
        )
    return AnnotationPanel(
        out_shards,
        names,
        is_binary=binary,
        annotation_metadata=metadata,
    )


def create_annotation(
    reference,
    annovec,
    annotation_name,
    is_binary=None,
    annotation_metadata=None,
) -> AnnotationPanel:
    name = str(annotation_name)
    if not name:
        raise ValueError("annotation_name must be non-empty")
    vec = _coerce_annovec(annovec, int(reference.num_snp))
    mat = sparse.csr_matrix(vec.reshape(-1, 1))
    if is_binary is not None:
        binary_arr = np.asarray(is_binary)
        if binary_arr.ndim != 0:
            raise ValueError("is_binary must be a scalar logical for create_annotation")
        binary = [bool(binary_arr)]
    else:
        binary = None
    metadata = "" if annotation_metadata is None else str(annotation_metadata)
    return create_annotations(
        reference,
        annotation_matrix=mat,
        annotation_names=[name],
        is_binary=binary,
        annotation_metadata=[metadata],
    )


def load_annotations(
    bed_paths,
    reference,
    annotation_metadata=None,
    annotation_metadata_paths=None,
) -> AnnotationPanel:
    if annotation_metadata is not None and annotation_metadata_paths is not None:
        raise ValueError("load_annotations accepts at most one of annotation_metadata and annotation_metadata_paths")

    if isinstance(bed_paths, (str, Path)):
        paths = [Path(bed_paths)]
    else:
        paths = [Path(p) for p in list(bed_paths)]
        if not paths:
            raise ValueError("bed_paths must be a non-empty list of BED files")
    for p in paths:
        if not p.is_file():
            raise FileNotFoundError(f"BED file not found: {p}")

    annomat, annonames, source_counts = _paint_binary_annotations(paths, reference)
    if annotation_metadata is not None:
        metadata = _coerce_string_vector(annotation_metadata, len(paths), "annotation_metadata")
    elif annotation_metadata_paths is not None:
        meta_paths = list(annotation_metadata_paths)
        if len(meta_paths) != len(paths):
            raise ValueError(
                f"annotation_metadata_paths length mismatch: expected {len(paths)}, got {len(meta_paths)}"
            )
        values = []
        for bed_path, meta_path, count in zip(paths, meta_paths, source_counts):
            if meta_path is None or str(meta_path) == "":
                values.append(_generated_metadata(bed_path, False, count))
            else:
                values.append(_read_sidecar_exact(meta_path))
        metadata = np.asarray(values, dtype=object)
    else:
        metadata = np.asarray(
            [_generated_metadata(p, False, count) for p, count in zip(paths, source_counts)],
            dtype=object,
        )

    return create_annotations(
        reference,
        annotation_matrix=annomat,
        annotation_names=annonames,
        is_binary=np.ones(len(annonames), dtype=bool),
        annotation_metadata=metadata,
    )


def load_annotation(
    path,
    reference,
    has_header: bool = False,
    value_columns=None,
    group_column=None,
    annotation_names=None,
    annotation_metadata=None,
    annotation_metadata_path=None,
) -> AnnotationPanel:
    if annotation_metadata is not None and annotation_metadata_path is not None:
        raise ValueError("load_annotation accepts at most one of annotation_metadata and annotation_metadata_path")
    if value_columns is not None and group_column is not None:
        raise ValueError("group_column and value_columns are mutually exclusive")
    if group_column is not None and annotation_names is not None:
        raise ValueError("annotation_names is invalid with group_column")
    if group_column is not None and annotation_metadata is not None:
        raise ValueError("annotation_metadata is invalid with group_column")
    if group_column is not None and annotation_metadata_path is not None:
        raise ValueError("annotation_metadata_path is invalid with group_column")

    source_path = Path(path)
    df, header_fields, row_base0 = _read_annotation_table(source_path, bool(has_header))
    num_columns = int(df.shape[1])
    if num_columns < 3:
        raise ValueError(f"{source_path}: annotation input must have at least 3 tab-separated columns")

    value_normalized = _normalize_value_columns(value_columns, header_fields, num_columns, source_path)
    value_columns0 = None if value_normalized is None else value_normalized[0]
    value_column_metadata = None if value_normalized is None else value_normalized[1]
    group_normalized = _normalize_group_column(group_column, header_fields, num_columns, source_path)
    group_column0 = None if group_normalized is None else group_normalized[0]
    group_column_metadata = None if group_normalized is None else group_normalized[1]
    chr_values, starts, ends = _validate_interval_columns(df, source_path, row_base0=row_base0)
    n = int(reference.num_snp)
    num_source_intervals = int(df.shape[0])

    if group_column0 is not None:
        group_values = df.iloc[:, group_column0].astype(str).to_numpy(dtype=object)
        empty = group_values == ""
        if empty.any():
            idx = int(np.flatnonzero(empty)[0])
            raise ValueError(f"{source_path}: row {row_base0 + idx + 1}: group_column values must be non-empty")
        names = [str(x) for x in pd.unique(group_values)]
        columns = []
        metadata_values = []
        # PERF: loop retained because grouped binary interval union is defined independently per group.
        for group_value in names:
            rows = group_values == group_value
            intervals_by_chr = binary_intervals_by_chr_from_arrays(chr_values[rows], starts[rows], ends[rows])
            columns.append(paint_binary_column(intervals_by_chr, reference))
            metadata_values.append(
                _generated_metadata(
                    source_path,
                    bool(has_header),
                    int(np.count_nonzero(rows)),
                    group_column=group_column_metadata,
                    group_value=group_value,
                )
            )
        annomat = sparse.hstack(columns, format="csr") if columns else sparse.csr_matrix((n, 0), dtype=np.float64)
        is_binary = np.ones(len(names), dtype=bool)
        metadata = np.asarray(metadata_values, dtype=object)
    elif value_columns0 is None and num_columns == 3:
        names = (
            _default_annotation_names(source_path, num_columns, header_fields, value_columns0)
            if annotation_names is None
            else _coerce_annonames(annotation_names)
        )
        if len(names) != 1:
            raise ValueError("annotation_names length mismatch: expected 1")
        intervals_by_chr = binary_intervals_by_chr_from_arrays(chr_values, starts, ends)
        annomat = paint_binary_column(intervals_by_chr, reference)
        is_binary = np.array([True], dtype=bool)
        if annotation_metadata_path is not None:
            metadata = np.asarray([_read_sidecar_exact(annotation_metadata_path)], dtype=object)
        elif annotation_metadata is not None:
            metadata = _coerce_string_vector(annotation_metadata, 1, "annotation_metadata")
        else:
            metadata = np.asarray(
                [_generated_metadata(source_path, bool(has_header), num_source_intervals)],
                dtype=object,
            )
    else:
        if value_columns0 is None:
            if num_columns == 4:
                value_columns0 = [3]
                value_column_metadata = [str(header_fields[3])] if header_fields is not None else [3]
            else:
                raise ValueError(f"{source_path}: input with five or more columns requires explicit value_columns or group_column")
        names = (
            _default_annotation_names(source_path, num_columns, header_fields, value_columns0)
            if annotation_names is None
            else _coerce_annonames(annotation_names)
        )
        if len(names) != len(value_columns0):
            raise ValueError(
                f"annotation_names length mismatch: expected {len(value_columns0)}, got {len(names)}"
            )

        validate_numeric_non_overlapping(
            chr_values,
            starts,
            ends,
            source_path,
            row_base0=row_base0,
        )
        values = _numeric_values(df, value_columns0, source_path, row_base0=row_base0)

        blocks = [
            sparse.csr_matrix((ref_shard.num_snp, len(value_columns0)), dtype=np.float64)
            for ref_shard in reference.shards
        ]
        source = pd.DataFrame({"chr": chr_values, "start": starts, "end": ends})
        for chr_label, grp in source.groupby("chr", sort=False):
            if chr_label in IGNORED_CHR:
                continue
            order = np.lexsort((grp["end"].to_numpy(), grp["start"].to_numpy()))
            intervals = np.column_stack(
                [
                    grp["start"].to_numpy(dtype=np.int64)[order],
                    grp["end"].to_numpy(dtype=np.int64)[order],
                ]
            )
            value_block = values[grp.index.to_numpy(dtype=np.int64)[order], :]
            for shard_idx, ref_shard in enumerate(reference.shards):
                if ref_shard.label != chr_label:
                    continue
                blocks[shard_idx] = paint_numeric_sparse(ref_shard.bp, intervals, value_block)
        annomat = sparse.vstack(blocks, format="csr")
        is_binary = np.zeros(len(value_columns0), dtype=bool)

        if annotation_metadata_path is not None:
            sidecar_lines = _read_column_metadata_sidecar(annotation_metadata_path, num_columns)
            metadata = np.asarray([sidecar_lines[i] for i in value_columns0], dtype=object)
        elif annotation_metadata is not None:
            metadata = _coerce_string_vector(
                annotation_metadata, len(value_columns0), "annotation_metadata"
            )
        else:
            metadata = np.asarray(
                [
                    _generated_metadata(
                        source_path,
                        bool(has_header),
                        num_source_intervals,
                        value_column=value_column,
                    )
                    for value_column in value_column_metadata
                ],
                dtype=object,
            )

    return create_annotations(
        reference,
        annotation_matrix=annomat,
        annotation_names=names,
        is_binary=is_binary,
        annotation_metadata=metadata,
    )


def save_annotations_cache(panel: AnnotationPanel, path) -> None:
    meta = {
        "schema": _CACHE_SCHEMA,
        "shard_labels": [s.label for s in panel.shards],
        "shard_checksums": [s.reference_checksum for s in panel.shards],
        "annonames": panel.annonames.tolist(),
        "is_binary": panel.is_binary.astype(bool).tolist(),
        "annotation_metadata": panel.annotation_metadata.tolist(),
    }
    arrays = {
        "_meta": np.frombuffer(json.dumps(meta, separators=(",", ":")).encode(), dtype=np.uint8),
    }

    for i, shard in enumerate(panel.shards):
        p = f"s{i}_"
        mat = shard.annomat.tocsr()
        arrays[p + "data"] = mat.data.astype(np.float64)
        arrays[p + "indices"] = mat.indices.astype(np.int32)
        arrays[p + "indptr"] = mat.indptr.astype(np.int32)
        arrays[p + "shape"] = np.asarray(mat.shape, dtype=np.int64)

    np.savez_compressed(path, **arrays)


def load_annotations_cache(path, shards=None) -> AnnotationPanel:
    with np.load(path, allow_pickle=False) as data:
        meta = json.loads(bytes(data["_meta"]).decode())
        schema = meta.get("schema")
        if schema not in {_CACHE_SCHEMA, _OLD_CACHE_SCHEMA}:
            raise ValueError(f"Unsupported annotations cache schema: {schema!r}")

        labels = list(meta.get("shard_labels", []))
        checksums = list(meta.get("shard_checksums", []))
        if len(labels) != len(checksums):
            raise ValueError(
                "Invalid annotations cache: shard_labels and shard_checksums length mismatch"
            )

        annonames = _coerce_annonames(meta.get("annonames", []))
        if schema == _OLD_CACHE_SCHEMA:
            is_binary = np.ones(len(annonames), dtype=bool)
            annotation_metadata = np.asarray([""] * len(annonames), dtype=object)
        else:
            is_binary = _coerce_is_binary_vector(meta.get("is_binary", []), len(annonames))
            annotation_metadata = _coerce_string_vector(
                meta.get("annotation_metadata", []),
                len(annonames),
                "annotation_metadata",
            )
        selected = validate_requested_shards(shards, labels, "load_annotations_cache")
        label_to_index = {label: i for i, label in enumerate(labels)}

        shard_objs = []
        for label in selected:
            i = label_to_index[label]
            p = f"s{i}_"
            mat = sparse.csr_matrix(
                (
                    np.asarray(data[p + "data"], dtype=np.float64),
                    np.asarray(data[p + "indices"], dtype=np.int32),
                    np.asarray(data[p + "indptr"], dtype=np.int32),
                ),
                shape=tuple(np.asarray(data[p + "shape"], dtype=np.int64).tolist()),
            )
            shard_objs.append(
                AnnotationShard._from_arrays(
                    label=label,
                    reference_checksum=checksums[i],
                    annomat=mat,
                )
            )

    return AnnotationPanel(
        shard_objs,
        annonames,
        is_binary=is_binary,
        annotation_metadata=annotation_metadata,
    )
