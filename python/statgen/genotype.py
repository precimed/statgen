from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ._bfile_utils import (
    FAM_PUBLIC_COLUMNS,
    _expected_bed_size,
    _parse_bim,
    _parse_fam,
    _parse_ploidy,
    _read_bed_rows_int8,
    _validate_bed,
)
from ._utils import allele_hash64, validate_requested_shards
from ._verbosity import report_variant_swaps
from ._id_match import match_id_panel, validate_matching_options
from ._variant_match import match_shard_numeric, numeric_variant_keys

_CACHE_SCHEMA = "genotype_cache/0.2"
_BED_MISSING_INT8 = np.int8(-1)
_SWAP_LUT = np.array([-1, 2, 1, 0], dtype=np.int8)
_HAPLOID_MODES = {"raw", "ploidy_scaled"}
_CACHE_ARRAY_FIELDS = (
    "is_present",
    "ploidy_male",
    "ploidy_female",
    "source_row0",
    "subject_present",
    "source_subject_row0",
    "fid",
    "iid",
    "father_id",
    "mother_id",
    "sex",
    "is_male",
    "is_female",
)


def _freeze(arr) -> np.ndarray:
    out = np.asarray(arr)
    out.setflags(write=False)
    return out


def _require_source_paths(prefix: str, label: str | None) -> tuple[Path, Path, Path, Path | None]:
    bed = Path(prefix + ".bed")
    bim = Path(prefix + ".bim")
    fam = Path(prefix + ".fam")
    for path in (bed, bim, fam):
        if not path.is_file():
            where = "" if label is None else f" for requested shard {label!r}"
            raise FileNotFoundError(f"Missing genotype source{where}: {path}")
    ploidy = Path(prefix + ".ploidy")
    return bed, bim, fam, ploidy if ploidy.is_file() else None


def _normalize_snp_indices(snp_indices, num_snp: int) -> np.ndarray:
    arr = np.asarray(snp_indices)
    if arr.ndim == 0:
        arr = arr.reshape(1)
    else:
        arr = arr.reshape(-1)
    if arr.size == 0:
        return np.array([], dtype=np.int64)
    if arr.dtype.kind not in {"i", "u"}:
        raise ValueError("GenotypePanel.fetch_genotypes: snp_indices must be integer indices")
    out = arr.astype(np.int64, copy=False)
    if np.any(out < 0) or np.any(out >= num_snp):
        bad = int(out[(out < 0) | (out >= num_snp)][0])
        raise ValueError(
            f"GenotypePanel.fetch_genotypes: SNP index {bad} is out of bounds for num_snp={num_snp}"
        )
    return out


def _normalize_haploid_mode(haploid_mode) -> str:
    mode = str(haploid_mode)
    if mode not in _HAPLOID_MODES:
        raise ValueError("GenotypePanel.fetch_genotypes: haploid_mode must be 'raw' or 'ploidy_scaled'")
    return mode

@dataclass(frozen=True)
class _SourceRecord:
    bed_path: Path
    bim_path: Path
    bed_file_size: int
    source_num_snp: int
    source_num_sample: int
    bim: pd.DataFrame
    fam: pd.DataFrame
    ploidy_male: np.ndarray
    ploidy_female: np.ndarray


def _load_source_record(prefix: str, label: str | None) -> _SourceRecord:
    bed, bim, fam, ploidy = _require_source_paths(prefix, label)
    source_bim = _parse_bim(bim)
    source_num_snp = int(source_bim.shape[0])
    source_bim = source_bim.copy()
    source_bim["line"] = np.arange(source_num_snp, dtype=np.int64) + 1
    source_bim["source_row0"] = np.arange(source_num_snp, dtype=np.int64)
    if label is not None:
        wrong = source_bim["chr"] != label
        if wrong.any():
            idx = int(wrong.idxmax())
            raise ValueError(
                f"{bim}:{int(source_bim.loc[idx, 'line'])}: sharded genotype BIM for {label!r} "
                f"contains chr {source_bim.loc[idx, 'chr']!r}"
            )
    source_bim["a1_hash64"] = allele_hash64(source_bim["a1"].to_numpy(dtype=object))
    source_bim["a2_hash64"] = allele_hash64(source_bim["a2"].to_numpy(dtype=object))
    source_fam = _parse_fam(fam)
    if label == "X" and ploidy is None:
        warnings.warn(
            f"{prefix}: chrX genotype source has no .ploidy sidecar; "
            "defaulting chrX rows to male/female ploidy (1, 2)",
            RuntimeWarning,
            stacklevel=2,
        )
    ploidy_male, ploidy_female = _parse_ploidy(
        ploidy,
        source_num_snp,
        source_bim["chr"].to_numpy(dtype=object),
    )
    bed_file_size = _validate_bed(bed, source_fam.shape[0], source_num_snp)
    return _SourceRecord(
        bed_path=bed,
        bim_path=bim,
        bed_file_size=bed_file_size,
        source_num_snp=source_num_snp,
        source_num_sample=int(source_fam.shape[0]),
        bim=source_bim,
        fam=source_fam,
        ploidy_male=ploidy_male,
        ploidy_female=ploidy_female,
    )


class GenotypeShard:
    def __init__(
        self,
        label: str,
        bed_path,
        bed_file_size,
        source_num_snp,
        source_num_sample,
        source_row0,
        subject_present,
        source_subject_row0,
        is_present,
        is_allele_swapped,
        ploidy_male,
        ploidy_female,
        reference_checksum,
    ):
        self._label = str(label)
        self._bed_path = Path(bed_path)
        self._bed_file_size = int(bed_file_size)
        self._source_num_snp = int(source_num_snp)
        self._source_num_sample = int(source_num_sample)
        if self._source_num_snp < 0 or self._source_num_sample < 0:
            raise ValueError("source_num_snp and source_num_sample must be non-negative")
        expected = _expected_bed_size(self._source_num_sample, self._source_num_snp)
        if self._bed_file_size != expected:
            raise ValueError(
                f"Invalid genotype metadata for shard {self._label}: bed_file_size "
                f"{self._bed_file_size} does not match expected {expected}"
            )

        self._source_row0 = np.asarray(source_row0, dtype=np.int64).reshape(-1)
        self._is_present = np.asarray(is_present, dtype=bool).reshape(-1)
        self._is_allele_swapped = np.asarray(is_allele_swapped, dtype=bool).reshape(-1)
        self._ploidy_male = np.asarray(ploidy_male, dtype=float).reshape(-1)
        self._ploidy_female = np.asarray(ploidy_female, dtype=float).reshape(-1)
        if not (
            self._source_row0.size
            == self._is_allele_swapped.size
            == self._is_present.size
            == self._ploidy_male.size
            == self._ploidy_female.size
        ):
            raise ValueError(f"GenotypeShard {self._label}: SNP-axis vector lengths must match")
        if not np.array_equal(self._is_present, self._source_row0 >= 0):
            raise ValueError(f"GenotypeShard {self._label}: is_present must match source_row0 >= 0")
        bad_source = (self._source_row0 < -1) | (self._source_row0 >= self._source_num_snp)
        if bad_source.any():
            raise ValueError(f"GenotypeShard {self._label}: source_row0 out of source BIM bounds")
        if np.any(self._is_allele_swapped & ~self._is_present):
            raise ValueError("is_allele_swapped requires is_present")
        absent = ~self._is_present
        if np.any(~np.isnan(self._ploidy_male[absent])) or np.any(~np.isnan(self._ploidy_female[absent])):
            raise ValueError(f"GenotypeShard {self._label}: absent SNPs must have NaN ploidy")
        present_ploidy = np.column_stack([self._ploidy_male[self._is_present], self._ploidy_female[self._is_present]])
        if present_ploidy.size and np.any(~np.isin(present_ploidy, [0.0, 1.0, 2.0])):
            raise ValueError(f"GenotypeShard {self._label}: present ploidy values must be 0, 1, or 2")

        self._subject_present = np.asarray(subject_present, dtype=bool).reshape(-1)
        self._source_subject_row0 = np.asarray(source_subject_row0, dtype=np.int64).reshape(-1)
        if self._subject_present.size != self._source_subject_row0.size:
            raise ValueError(f"GenotypeShard {self._label}: sample-axis vector lengths must match")
        if not np.array_equal(self._subject_present, self._source_subject_row0 >= 0):
            raise ValueError(f"GenotypeShard {self._label}: subject_present must match source_subject_row0 >= 0")
        bad_subject = (self._source_subject_row0 < -1) | (
            self._source_subject_row0 >= self._source_num_sample
        )
        if bad_subject.any():
            raise ValueError(f"GenotypeShard {self._label}: source_subject_row0 out of source FAM bounds")

        self._reference_checksum = str(reference_checksum)
        for name in (
            "_source_row0",
            "_is_present",
            "_is_allele_swapped",
            "_ploidy_male",
            "_ploidy_female",
            "_subject_present",
            "_source_subject_row0",
        ):
            getattr(self, name).setflags(write=False)

    @property
    def label(self) -> str:
        return self._label

    @property
    def chr(self) -> str:
        return self._label

    @property
    def num_snp(self) -> int:
        return int(self._is_present.size)

    @property
    def bed_path(self) -> Path:
        return self._bed_path

    @property
    def bed_file_size(self) -> int:
        return self._bed_file_size

    @property
    def source_num_snp(self) -> int:
        return self._source_num_snp

    @property
    def source_num_sample(self) -> int:
        return self._source_num_sample

    @property
    def source_row0(self) -> np.ndarray:
        return self._source_row0

    @property
    def subject_present(self) -> np.ndarray:
        return self._subject_present

    @property
    def source_subject_row0(self) -> np.ndarray:
        return self._source_subject_row0

    @property
    def is_present(self) -> np.ndarray:
        return self._is_present

    @property
    def is_allele_swapped(self) -> np.ndarray:
        return self._is_allele_swapped

    @property
    def ploidy_male(self) -> np.ndarray:
        return self._ploidy_male

    @property
    def ploidy_female(self) -> np.ndarray:
        return self._ploidy_female

    @property
    def reference_checksum(self) -> str:
        return self._reference_checksum


class GenotypePanel:
    def __init__(
        self,
        shards: list[GenotypeShard],
        *,
        fid,
        iid,
        father_id,
        mother_id,
        sex,
        source_layout: str,
    ):
        if source_layout not in {"non_sharded", "sharded"}:
            raise ValueError("source_layout must be 'non_sharded' or 'sharded'")
        self._source_layout = source_layout
        self._shards = list(shards)
        if not self._shards:
            raise ValueError("GenotypePanel requires at least one shard")
        self._fid = _freeze(np.asarray(fid, dtype=object).reshape(-1))
        self._iid = _freeze(np.asarray(iid, dtype=object).reshape(-1))
        self._father_id = _freeze(np.asarray(father_id, dtype=object).reshape(-1))
        self._mother_id = _freeze(np.asarray(mother_id, dtype=object).reshape(-1))
        self._sex = _freeze(np.asarray(sex, dtype=np.int8).reshape(-1))
        if not (
            self._fid.size
            == self._iid.size
            == self._father_id.size
            == self._mother_id.size
            == self._sex.size
        ):
            raise ValueError("GenotypePanel FAM vector lengths must match")
        if np.any(~np.isin(self._sex, [0, 1, 2])):
            raise ValueError("GenotypePanel sex values must be 0, 1, or 2")

        self._shard_offsets = []
        pos = 0
        for shard in self._shards:
            if shard.subject_present.size != self._fid.size:
                raise ValueError(
                    f"Genotype shard {shard.label}: subject_present length does not match panel sample axis"
                )
            self._shard_offsets.append(
                {"shard_label": shard.label, "start0": pos, "stop0": pos + shard.num_snp}
            )
            pos += shard.num_snp
        self._num_snp = int(pos)

    @property
    def shards(self) -> list[GenotypeShard]:
        return list(self._shards)

    @property
    def shard_offsets(self) -> list[dict]:
        return list(self._shard_offsets)

    @property
    def source_layout(self) -> str:
        return self._source_layout

    @property
    def num_snp(self) -> int:
        return self._num_snp

    @property
    def is_present(self) -> np.ndarray:
        if not self._shards:
            return np.array([], dtype=bool)
        return np.concatenate([s.is_present for s in self._shards])

    @property
    def is_allele_swapped(self) -> np.ndarray:
        if not self._shards:
            return np.array([], dtype=bool)
        return np.concatenate([s.is_allele_swapped for s in self._shards])

    @property
    def ploidy_male(self) -> np.ndarray:
        if not self._shards:
            return np.array([], dtype=float)
        return np.concatenate([s.ploidy_male for s in self._shards])

    @property
    def ploidy_female(self) -> np.ndarray:
        if not self._shards:
            return np.array([], dtype=float)
        return np.concatenate([s.ploidy_female for s in self._shards])

    @property
    def source_row0(self) -> np.ndarray:
        if not self._shards:
            return np.array([], dtype=np.int64)
        return np.concatenate([s.source_row0 for s in self._shards])

    @property
    def num_sample(self) -> int:
        return int(self._fid.size)

    @property
    def fid(self) -> np.ndarray:
        return self._fid

    @property
    def iid(self) -> np.ndarray:
        return self._iid

    @property
    def father_id(self) -> np.ndarray:
        return self._father_id

    @property
    def mother_id(self) -> np.ndarray:
        return self._mother_id

    @property
    def sex(self) -> np.ndarray:
        return self._sex

    @property
    def is_male(self) -> np.ndarray:
        return self._sex == 1

    @property
    def is_female(self) -> np.ndarray:
        return self._sex == 2

    def is_subject_present(self, shard) -> np.ndarray:
        label = str(shard)
        for s in self._shards:
            if s.label == label:
                return s.subject_present
        raise ValueError(f"GenotypePanel.is_subject_present: requested shard {label!r} is not loaded")

    def select_shards(self, shards) -> "GenotypePanel":
        available = [s.label for s in self._shards]
        selected = validate_requested_shards(shards, available, "GenotypePanel.select_shards")
        by_label = {s.label: s for s in self._shards}
        return GenotypePanel(
            [by_label[label] for label in selected],
            fid=self._fid,
            iid=self._iid,
            father_id=self._father_id,
            mother_id=self._mother_id,
            sex=self._sex,
            source_layout=self._source_layout,
        )

    def _addressed_shards(self, snp_indices: np.ndarray) -> list[tuple[GenotypeShard, np.ndarray, np.ndarray]]:
        addressed = []
        for shard, offset in zip(self._shards, self._shard_offsets):
            start = int(offset["start0"])
            stop = int(offset["stop0"])
            cols = np.flatnonzero((snp_indices >= start) & (snp_indices < stop))
            if cols.size:
                addressed.append((shard, cols, snp_indices[cols] - start))
        return addressed

    def _validate_requested_present(self, addressed) -> None:
        start_by_label = {str(o["shard_label"]): int(o["start0"]) for o in self._shard_offsets}
        for shard, _cols, local_indices in addressed:
            present = shard.is_present[local_indices]
            if not np.all(present):
                j = int(np.flatnonzero(~present)[0])
                panel_index = start_by_label[shard.label] + int(local_indices[j])
                raise ValueError(
                    "GenotypePanel.fetch_genotypes: requested SNP "
                    f"{panel_index} in shard {shard.label!r} is not present in the genotype source"
                )

    def fetch_genotypes_int8(self, snp_indices):
        indices = _normalize_snp_indices(snp_indices, self._num_snp)
        out = np.full((self.num_sample, indices.size), _BED_MISSING_INT8, dtype=np.int8)
        if indices.size == 0:
            return out

        addressed = self._addressed_shards(indices)
        self._validate_requested_present(addressed)

        for shard, cols, local_indices in addressed:
            effective_bed = shard.bed_path
            _validate_bed(effective_bed, shard.source_num_sample, shard.source_num_snp)
            source_rows = shard.source_row0[local_indices]
            source_geno = _read_bed_rows_int8(
                effective_bed,
                source_rows,
                shard.source_num_sample,
            )
            sc = np.flatnonzero(shard.is_allele_swapped[local_indices])
            if sc.size:
                source_geno[:, sc] = _SWAP_LUT[source_geno[:, sc] + 1]
            panel_rows = np.flatnonzero(shard.subject_present)
            if panel_rows.size:
                source_subject_rows = shard.source_subject_row0[panel_rows]
                out[np.ix_(panel_rows, cols)] = source_geno[source_subject_rows, :]
        return out

    def fetch_genotypes(self, snp_indices, haploid_mode="raw"):
        mode = _normalize_haploid_mode(haploid_mode)
        geno_int8 = self.fetch_genotypes_int8(snp_indices)
        geno = geno_int8.astype(np.float64)
        geno[geno_int8 == _BED_MISSING_INT8] = np.nan
        if mode == "ploidy_scaled":
            self._apply_ploidy_scaled(geno, _normalize_snp_indices(snp_indices, self._num_snp))
        return geno

    def _apply_ploidy_scaled(self, geno: np.ndarray, indices: np.ndarray) -> None:
        if indices.size == 0:
            return
        is_male = self._sex == 1
        is_female = self._sex == 2
        is_unknown = self._sex == 0
        addressed = self._addressed_shards(indices)
        for shard, cols, local_indices in addressed:
            ploidy_male = shard.ploidy_male[local_indices]
            ploidy_female = shard.ploidy_female[local_indices]
            differs = ploidy_male != ploidy_female
            unknown_present = is_unknown & shard.subject_present
            if (
                np.any(unknown_present)
                and np.any(differs)
                and np.any(np.isfinite(geno[np.ix_(unknown_present, cols[differs])]))
            ):
                raise ValueError(
                    "GenotypePanel.fetch_genotypes: haploid_mode='ploidy_scaled' "
                    "requires known FAM sex when male and female ploidy differ"
                )
            if np.any(is_male):
                geno[np.ix_(is_male, cols)] *= ploidy_male.reshape(1, -1) / 2.0
            if np.any(is_female):
                geno[np.ix_(is_female, cols)] *= ploidy_female.reshape(1, -1) / 2.0
            if np.any(is_unknown):
                same_ploidy = ploidy_male.reshape(1, -1)
                geno[np.ix_(is_unknown, cols)] *= same_ploidy / 2.0

    def save_cache(self, path) -> None:
        save_genotype_cache(self, path)


def _fam_equal(lhs: pd.DataFrame, rhs: pd.DataFrame) -> bool:
    return lhs[FAM_PUBLIC_COLUMNS].reset_index(drop=True).equals(
        rhs[FAM_PUBLIC_COLUMNS].reset_index(drop=True)
    )


def _all_subjects_present(source_num_sample: int) -> tuple[np.ndarray, np.ndarray]:
    return (
        np.ones(source_num_sample, dtype=bool),
        np.arange(source_num_sample, dtype=np.int64),
    )


def _map_chrx_subjects(chrx_fam: pd.DataFrame, panel_fam: pd.DataFrame, where: str) -> tuple[np.ndarray, np.ndarray]:
    num_sample = int(panel_fam.shape[0])
    if _fam_equal(chrx_fam, panel_fam):
        return _all_subjects_present(num_sample)

    panel = panel_fam.copy()
    panel["panel_row0"] = np.arange(num_sample, dtype=np.int64)
    src = chrx_fam.copy()
    src["source_row0"] = np.arange(src.shape[0], dtype=np.int64)
    merged = src.merge(panel, on=["fid", "iid"], how="left", sort=False, suffixes=("_src", "_panel"))
    missing = merged["panel_row0"].isna()
    if missing.any():
        idx = int(missing.idxmax())
        raise ValueError(
            f"{where}: chrX FAM subject ({merged.loc[idx, 'fid']!r}, {merged.loc[idx, 'iid']!r}) "
            "is absent from the autosomal sample axis"
        )

    for col in ("father_id", "mother_id", "sex"):
        bad = merged[f"{col}_src"] != merged[f"{col}_panel"]
        if bad.any():
            idx = int(bad.idxmax())
            raise ValueError(
                f"{where}: chrX FAM metadata mismatch for subject "
                f"({merged.loc[idx, 'fid']!r}, {merged.loc[idx, 'iid']!r})"
            )

    panel_rows = merged["panel_row0"].to_numpy(dtype=np.int64)
    source_rows = merged["source_row0"].to_numpy(dtype=np.int64)
    subject_present = np.zeros(num_sample, dtype=bool)
    source_subject_row0 = np.full(num_sample, -1, dtype=np.int64)
    subject_present[panel_rows] = True
    source_subject_row0[panel_rows] = source_rows
    return subject_present, source_subject_row0


def _build_shard(ref_shard, source: _SourceRecord, subject_present, source_subject_row0, allow_indel_swaps, id_match=None) -> GenotypeShard:
    if id_match is not None:
        source_bim = source.bim
        local_match, is_allele_swapped = id_match
    elif (source_bim := source.bim.loc[source.bim["chr"] == ref_shard.label]).empty:
        local_match = np.full(int(ref_shard.num_snp), -1, dtype=np.int64)
        is_allele_swapped = np.zeros(int(ref_shard.num_snp), dtype=bool)
    else:
        ref_keys = numeric_variant_keys(ref_shard.bp, ref_shard.a1_hash64, ref_shard.a2_hash64)
        src_keys = numeric_variant_keys(
            source_bim["bp"].to_numpy(dtype=np.int64),
            source_bim["a1_hash64"].to_numpy(dtype=np.uint64),
            source_bim["a2_hash64"].to_numpy(dtype=np.uint64),
        )
        match = match_shard_numeric(
            ref_keys, src_keys, ref_shard.label, source_name="genotype",
            allow_indel_swaps=allow_indel_swaps,
        )
        local_match = match.source_rows
        is_allele_swapped = match.is_allele_swapped
        report_variant_swaps(source.bim_path, ref_shard.label, np.count_nonzero(is_allele_swapped),
                             match.num_strand_ambiguous, match.num_disabled_indels)

    is_present = local_match >= 0
    source_row0 = np.full(int(ref_shard.num_snp), -1, dtype=np.int64)
    ploidy_male = np.full(int(ref_shard.num_snp), np.nan, dtype=float)
    ploidy_female = np.full(int(ref_shard.num_snp), np.nan, dtype=float)
    if is_present.any():
        source_rows = source_bim["source_row0"].to_numpy(dtype=np.int64)[local_match[is_present]]
        source_row0[is_present] = source_rows
        ploidy_male[is_present] = source.ploidy_male[source_rows]
        ploidy_female[is_present] = source.ploidy_female[source_rows]

    return GenotypeShard(
        label=ref_shard.label,
        bed_path=source.bed_path,
        bed_file_size=source.bed_file_size,
        source_num_snp=source.source_num_snp,
        source_num_sample=source.source_num_sample,
        source_row0=source_row0,
        subject_present=subject_present,
        source_subject_row0=source_subject_row0,
        is_present=is_present,
        is_allele_swapped=is_allele_swapped,
        ploidy_male=ploidy_male,
        ploidy_female=ploidy_female,
        reference_checksum=ref_shard.checksum,
    )


def _validate_reference_for_genotype(reference) -> list:
    ref_shards = list(getattr(reference, "shards", []))
    if not ref_shards:
        raise ValueError("load_genotype requires a reference with at least one shard")
    labels = [s.label for s in ref_shards]
    validate_requested_shards(labels, labels, "load_genotype reference")
    return ref_shards


def load_genotype(bfile_prefix, reference, allow_indel_swaps=False, match_by="genomic_position") -> GenotypePanel:
    allow_indel_swaps = validate_matching_options(match_by, allow_indel_swaps)
    ref_shards = _validate_reference_for_genotype(reference)
    prefix = str(bfile_prefix)
    ref_labels = [s.label for s in ref_shards]

    if "@" in prefix:
        source_layout = "sharded"
        sources = {
            label: _load_source_record(prefix.replace("@", label), label)
            for label in ref_labels
        }
    else:
        source_layout = "non_sharded"
        shared = _load_source_record(prefix, None)
        if (
            "X" in ref_labels
            and (shared.bim["chr"] == "X").any()
            and Path(prefix + ".ploidy").is_file() is False
        ):
            warnings.warn(
                f"{prefix}: chrX genotype source has no .ploidy sidecar; "
                "defaulting chrX rows to male/female ploidy (1, 2)",
                RuntimeWarning,
                stacklevel=2,
            )
        sources = {label: shared for label in ref_labels}

    subject_maps: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    if source_layout == "non_sharded":
        panel_fam = next(iter(sources.values())).fam
        for label in ref_labels:
            if sources[label].bim.loc[sources[label].bim["chr"] == label].empty:
                warnings.warn(
                    f"{prefix}: no source BIM rows for requested reference shard {label!r}; marking shard absent",
                    RuntimeWarning,
                    stacklevel=2,
                )
            subject_maps[label] = _all_subjects_present(sources[label].source_num_sample)
    else:
        autosomal_labels = [label for label in ref_labels if label != "X"]
        if autosomal_labels:
            panel_fam = sources[autosomal_labels[0]].fam
            for label in autosomal_labels[1:]:
                if not _fam_equal(sources[label].fam, panel_fam):
                    raise ValueError(f"load_genotype: autosomal FAM mismatch for shard {label!r}")
            for label in autosomal_labels:
                subject_maps[label] = _all_subjects_present(sources[label].source_num_sample)
            if "X" in sources:
                subject_maps["X"] = _map_chrx_subjects(sources["X"].fam, panel_fam, "load_genotype")
        else:
            panel_fam = sources["X"].fam
            subject_maps["X"] = _all_subjects_present(sources["X"].source_num_sample)

    id_matches = {}
    if match_by == "id":
        records = [shared] if source_layout == "non_sharded" else list(sources.values())
        inputs = [
            (record.bim_path, record.bim["snp"].to_numpy(dtype=object),
             record.bim["a1_hash64"].to_numpy(dtype=np.uint64),
             record.bim["a2_hash64"].to_numpy(dtype=np.uint64),
             record.bim["chr"].to_numpy(dtype=object))
            for record in records
        ]
        rows, swapped = match_id_panel(reference, inputs, chromosome_scope=True)
        source_offsets = {}
        offset = 0
        for record in records:
            source_offsets[id(record)] = offset
            offset += len(record.bim)
        for shard, off in zip(ref_shards, reference.shard_offsets):
            sl = slice(off["start0"], off["stop0"])
            local = rows[sl].copy()
            local[local >= 0] -= source_offsets[id(sources[shard.label])]
            id_matches[shard.label] = (local, swapped[sl])

    shards = [
        _build_shard(ref_shard, sources[ref_shard.label], *subject_maps[ref_shard.label], allow_indel_swaps,
                     id_matches.get(ref_shard.label))
        for ref_shard in ref_shards
    ]
    return GenotypePanel(
        shards,
        fid=panel_fam["fid"].to_numpy(dtype=object),
        iid=panel_fam["iid"].to_numpy(dtype=object),
        father_id=panel_fam["father_id"].to_numpy(dtype=object),
        mother_id=panel_fam["mother_id"].to_numpy(dtype=object),
        sex=panel_fam["sex"].to_numpy(dtype=np.int8),
        source_layout=source_layout,
    )


def save_genotype_cache(panel: GenotypePanel, path) -> None:
    meta = {
        "schema": _CACHE_SCHEMA,
        "n_shards": len(panel.shards),
        "shard_labels": [s.label for s in panel.shards],
        "shard_checksums": [s.reference_checksum for s in panel.shards],
        "shard_start0": [int(o["start0"]) for o in panel.shard_offsets],
        "shard_stop0": [int(o["stop0"]) for o in panel.shard_offsets],
        "source_layout": panel.source_layout,
        "bed_paths": [str(s.bed_path) for s in panel.shards],
        "bed_file_sizes": [s.bed_file_size for s in panel.shards],
        "source_num_snp": [s.source_num_snp for s in panel.shards],
        "source_num_sample": [s.source_num_sample for s in panel.shards],
        "num_sample": panel.num_sample,
    }
    arrays = {
        "_meta": np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8),
        "is_present": panel.is_present.astype(bool, copy=False),
        "is_allele_swapped": panel.is_allele_swapped,
        "ploidy_male": panel.ploidy_male.astype(float, copy=False),
        "ploidy_female": panel.ploidy_female.astype(float, copy=False),
        "source_row0": panel.source_row0.astype(np.int64, copy=False),
        "subject_present": np.column_stack([s.subject_present for s in panel.shards]).astype(bool, copy=False),
        "source_subject_row0": np.column_stack([s.source_subject_row0 for s in panel.shards]).astype(np.int64, copy=False),
        "fid": np.asarray(panel.fid, dtype=str),
        "iid": np.asarray(panel.iid, dtype=str),
        "father_id": np.asarray(panel.father_id, dtype=str),
        "mother_id": np.asarray(panel.mother_id, dtype=str),
        "sex": panel.sex.astype(np.int8, copy=False),
        "is_male": panel.is_male.astype(bool, copy=False),
        "is_female": panel.is_female.astype(bool, copy=False),
    }
    np.savez_compressed(path, **arrays)


def _load_cache_meta(data) -> dict:
    if "_meta" not in data.files:
        raise ValueError("Invalid genotype cache: missing _meta")
    try:
        meta = json.loads(bytes(data["_meta"]).decode())
    except (UnicodeDecodeError, json.JSONDecodeError, TypeError) as exc:
        raise ValueError("Invalid genotype cache: malformed _meta") from exc
    if not isinstance(meta, dict):
        raise ValueError("Invalid genotype cache: _meta must be a JSON object")
    return meta


def _validate_cache_meta(meta: dict, data) -> None:
    schema = meta.get("schema")
    if schema not in {_CACHE_SCHEMA, "genotype_cache/0.1"}:
        raise ValueError(f"Unsupported genotype cache schema: {schema!r}")
    required = _CACHE_ARRAY_FIELDS + (("is_allele_swapped",) if schema == _CACHE_SCHEMA else ())
    missing = [name for name in required if name not in data.files]
    if missing:
        raise ValueError(f"Invalid genotype cache: missing array field {missing[0]!r}")
    n_shards = int(meta.get("n_shards", -1))
    if n_shards <= 0:
        raise ValueError("Invalid genotype cache: n_shards must be positive")
    labels = list(meta.get("shard_labels", []))
    checksums = list(meta.get("shard_checksums", []))
    starts = list(meta.get("shard_start0", []))
    stops = list(meta.get("shard_stop0", []))
    vector_fields = (
        checksums,
        starts,
        stops,
        list(meta.get("bed_paths", [])),
        list(meta.get("bed_file_sizes", [])),
        list(meta.get("source_num_snp", [])),
        list(meta.get("source_num_sample", [])),
    )
    if len(labels) != n_shards or any(len(v) != n_shards for v in vector_fields):
        raise ValueError("Invalid genotype cache: shard metadata length mismatch")
    validate_requested_shards(labels, labels, "load_genotype_cache")
    source_layout = meta.get("source_layout")
    if source_layout not in {"non_sharded", "sharded"}:
        raise ValueError("Invalid genotype cache: source_layout must be 'non_sharded' or 'sharded'")

    is_present = data["is_present"]
    if schema == _CACHE_SCHEMA:
        swapped = data["is_allele_swapped"]
        if swapped.dtype != np.dtype(bool) or swapped.shape != is_present.shape:
            raise ValueError("Invalid genotype cache: is_allele_swapped must be a logical SNP-axis vector")
        if np.any(swapped & ~is_present.astype(bool)):
            raise ValueError("Invalid genotype cache: is_allele_swapped requires is_present")
    source_row0 = data["source_row0"]
    if is_present.ndim != 1 or source_row0.ndim != 1 or is_present.shape[0] != source_row0.shape[0]:
        raise ValueError("Invalid genotype cache: SNP-axis vector lengths mismatch")
    if data["ploidy_male"].shape != is_present.shape or data["ploidy_female"].shape != is_present.shape:
        raise ValueError("Invalid genotype cache: ploidy vector lengths mismatch")
    if starts and (starts[0] != 0 or stops[-1] != is_present.shape[0]):
        raise ValueError("Invalid genotype cache: shard offsets do not cover SNP-axis vectors")
    if any(int(starts[i]) >= int(stops[i]) for i in range(n_shards)):
        raise ValueError("Invalid genotype cache: shard offsets must be non-empty half-open ranges")
    if any(int(starts[i]) != int(stops[i - 1]) for i in range(1, n_shards)):
        raise ValueError("Invalid genotype cache: shard offsets must be contiguous")

    num_sample = int(meta.get("num_sample", -1))
    if num_sample < 0:
        raise ValueError("Invalid genotype cache: num_sample must be non-negative")
    for name in ("fid", "iid", "father_id", "mother_id", "sex", "is_male", "is_female"):
        if data[name].shape != (num_sample,):
            raise ValueError(f"Invalid genotype cache: {name} length mismatch")
    sex = data["sex"]
    if np.any(~np.isin(sex, [0, 1, 2])):
        raise ValueError("Invalid genotype cache: sex values must be 0, 1, or 2")
    if not np.array_equal(data["is_male"], sex == 1):
        raise ValueError("Invalid genotype cache: is_male does not match sex")
    if not np.array_equal(data["is_female"], sex == 2):
        raise ValueError("Invalid genotype cache: is_female does not match sex")
    if data["subject_present"].shape != (num_sample, n_shards):
        raise ValueError("Invalid genotype cache: subject_present shape mismatch")
    if data["source_subject_row0"].shape != (num_sample, n_shards):
        raise ValueError("Invalid genotype cache: source_subject_row0 shape mismatch")

    for i in range(n_shards):
        expected = _expected_bed_size(meta["source_num_sample"][i], meta["source_num_snp"][i])
        if int(meta["bed_file_sizes"][i]) != expected:
            raise ValueError(f"Invalid genotype cache: bed_file_size mismatch for shard {labels[i]!r}")


def _resolve_bed_paths(source_layout, labels, stored_paths, source_num_snp, source_num_sample, bed_path):
    if bed_path is None:
        return {label: Path(path) for label, path in zip(labels, stored_paths)}

    override = str(bed_path)
    if "@" in override:
        if source_layout == "non_sharded":
            raise ValueError(
                "load_genotype_cache: @ override incompatible with non-sharded panel metadata"
            )
        return {label: Path(override.replace("@", label)) for label in labels}

    if source_layout == "sharded" and labels:
        if len(set(source_num_snp)) != 1 or len(set(source_num_sample)) != 1:
            raise ValueError(
                "load_genotype_cache: flat bed_path override requires equal "
                "source_num_snp and source_num_sample across selected shards"
            )
    flat = Path(override)
    return {label: flat for label in labels}


def load_genotype_cache(path, shards=None, bed_path=None) -> GenotypePanel:
    """Load cached metadata, optionally replacing BED paths for selected shards.

    Missing or incorrectly sized stored BED files emit a warning; an invalid
    override raises an error. BED headers are checked when genotypes are read.
    """
    with np.load(path, allow_pickle=False) as data:
        meta = _load_cache_meta(data)
        _validate_cache_meta(meta, data)
        labels = list(meta["shard_labels"])
        selected = validate_requested_shards(shards, labels, "load_genotype_cache")
        label_to_index = {label: i for i, label in enumerate(labels)}

        indices = [label_to_index[label] for label in selected]
        selected_paths = _resolve_bed_paths(
            meta["source_layout"], selected, [meta["bed_paths"][i] for i in indices],
            [meta["source_num_snp"][i] for i in indices],
            [meta["source_num_sample"][i] for i in indices], bed_path,
        )
        bad = []
        for label, i in zip(selected, indices):
            bed = selected_paths[label]
            try:
                if not bed.is_file() or bed.stat().st_size != int(meta["bed_file_sizes"][i]):
                    bad.append(str(bed))
            except OSError:
                bad.append(str(bed))
        if bad:
            message = "Missing or incorrectly sized BED files: " + ", ".join(dict.fromkeys(bad))
            if bed_path is not None:
                raise ValueError(message)
            warnings.warn(message, UserWarning, stacklevel=2)
        swap_flags = (data["is_allele_swapped"] if meta["schema"] == _CACHE_SCHEMA
                      else np.zeros(data["is_present"].size, dtype=bool))
        shard_objs = []
        for label in selected:
            i = label_to_index[label]
            start = int(meta["shard_start0"][i])
            stop = int(meta["shard_stop0"][i])
            shard_objs.append(
                GenotypeShard(
                    label=label,
                    bed_path=selected_paths[label],
                    bed_file_size=int(meta["bed_file_sizes"][i]),
                    source_num_snp=int(meta["source_num_snp"][i]),
                    source_num_sample=int(meta["source_num_sample"][i]),
                    source_row0=data["source_row0"][start:stop],
                    subject_present=data["subject_present"][:, i],
                    source_subject_row0=data["source_subject_row0"][:, i],
                    is_present=data["is_present"][start:stop],
                    is_allele_swapped=swap_flags[start:stop],
                    ploidy_male=data["ploidy_male"][start:stop],
                    ploidy_female=data["ploidy_female"][start:stop],
                    reference_checksum=meta["shard_checksums"][i],
                )
            )

        return GenotypePanel(
            shard_objs,
            fid=data["fid"],
            iid=data["iid"],
            father_id=data["father_id"],
            mother_id=data["mother_id"],
            sex=data["sex"],
            source_layout=meta["source_layout"],
        )
