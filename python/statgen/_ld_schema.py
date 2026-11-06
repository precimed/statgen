from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np


MANIFEST_SCHEMA = "1.0"
SHARD_SCHEMA = "1.0"
PY_RUNTIME_FORMAT = "python_npz_csc32"
MAT_RUNTIME_FORMAT = "matlab_mat_sparse_double"
NPZ_FORMAT = "statgen_ld_npz_csc32"
VALID_CHRX_SEX = {"female", "male", "combined"}
NPZ_REQUIRED = {"data", "indices", "indptr", "shape", "a1freq", "metadata"}


def read_manifest(path: Path, expected_runtime: str | None) -> dict:
    require_file(path)
    with open(path, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    if manifest.get("object_type") != "ld_panel_manifest":
        raise ValueError(f"{path}: object_type must be 'ld_panel_manifest'")
    if manifest.get("schema_version") != MANIFEST_SCHEMA:
        raise ValueError(f"{path}: unsupported LD manifest schema_version {manifest.get('schema_version')!r}")
    runtime = manifest.get("runtime_format")
    if runtime not in {PY_RUNTIME_FORMAT, MAT_RUNTIME_FORMAT}:
        raise ValueError(f"{path}: unsupported LD runtime_format {runtime!r}")
    if expected_runtime is not None and runtime != expected_runtime:
        raise ValueError(f"{path}: expected runtime_format {expected_runtime!r}, got {runtime!r}")
    validate_plain_relative_filename(
        manifest.get("reference_cache"),
        f"{path}: reference_cache",
    )
    reference_cache_md5 = manifest.get("reference_cache_md5")
    if not isinstance(reference_cache_md5, str) or len(reference_cache_md5) != 32:
        raise ValueError(f"{path}: reference_cache_md5 must be a lowercase MD5 hex string")
    shards = manifest.get("shards")
    if not isinstance(shards, list) or not shards:
        raise ValueError(f"{path}: shards must be a non-empty list")

    seen = set()
    reference_bim_by_chr = {}
    for i, entry in enumerate(shards):
        validate_manifest_entry(entry, f"{path}:shards[{i}]")
        key = (entry["chr"], entry["sex"])
        if key in seen:
            raise ValueError(f"{path}: duplicate manifest entry for chr {key[0]} sex {key[1]!r}")
        seen.add(key)
        prior_reference_bim = reference_bim_by_chr.setdefault(entry["chr"], entry["reference_bim"])
        if entry["reference_bim"] != prior_reference_bim:
            raise ValueError(f"{path}: manifest entries for chr {entry['chr']} must share reference_bim")
    return manifest


def validate_manifest_entry(entry: dict, where: str) -> None:
    required = {
        "chr",
        "sex",
        "file",
        "file_md5",
        "num_snp",
        "nnz",
        "reference_checksum",
        "reference_bim",
    }
    missing = sorted(required.difference(entry))
    if missing:
        raise ValueError(f"{where}: missing required fields: {', '.join(missing)}")
    validate_chr_sex(entry["chr"], entry["sex"], where)
    if not isinstance(entry["file"], str) or entry["file"] == "":
        raise ValueError(f"{where}: file must be a non-empty relative path")
    if Path(entry["file"]).is_absolute() or ".." in Path(entry["file"]).parts:
        raise ValueError(f"{where}: file must be a relative path inside the LD distribution")
    validate_positive_int(entry["num_snp"], f"{where}: num_snp")
    validate_nonnegative_int(entry["nnz"], f"{where}: nnz")
    if not isinstance(entry["reference_checksum"], str) or entry["reference_checksum"] == "":
        raise ValueError(f"{where}: reference_checksum must be a non-empty string")
    if not isinstance(entry["file_md5"], str) or len(entry["file_md5"]) != 32:
        raise ValueError(f"{where}: file_md5 must be a lowercase MD5 hex string")
    validate_plain_relative_filename(entry["reference_bim"], f"{where}: reference_bim")


def validate_shard_metadata(meta: dict, path: Path, expected_format: str) -> None:
    required = {
        "object_type",
        "schema_version",
        "format",
        "chr",
        "sex",
        "num_snp",
        "nnz",
        "matrix",
        "diagonal",
        "value",
        "reference_checksum",
        "reference_bim",
        "num_monomorphic_snps",
    }
    missing = sorted(required.difference(meta))
    if missing:
        raise ValueError(f"{path}: metadata missing required fields: {', '.join(missing)}")
    if meta["object_type"] != "ld_shard":
        raise ValueError(f"{path}: metadata object_type must be 'ld_shard'")
    if meta["schema_version"] != SHARD_SCHEMA:
        raise ValueError(f"{path}: unsupported metadata schema_version {meta['schema_version']!r}")
    if meta["format"] != expected_format:
        raise ValueError(f"{path}: metadata format must be {expected_format!r}")
    validate_chr_sex(meta["chr"], meta["sex"], f"{path}: metadata")
    validate_positive_int(meta["num_snp"], f"{path}: metadata num_snp")
    validate_nonnegative_int(meta["nnz"], f"{path}: metadata nnz")
    validate_nonnegative_int(
        meta["num_monomorphic_snps"],
        f"{path}: metadata num_monomorphic_snps",
    )
    if int(meta["num_monomorphic_snps"]) > int(meta["num_snp"]):
        raise ValueError(f"{path}: metadata num_monomorphic_snps must not exceed num_snp")
    if int(meta["nnz"]) >= 2**31:
        raise ValueError(f"{path}: CSC32 metadata nnz must be < 2^31")
    if meta["matrix"] != "symmetric":
        raise ValueError(f"{path}: metadata matrix must be 'symmetric'")
    if meta["diagonal"] != "explicit_unit":
        raise ValueError(f"{path}: metadata diagonal must be 'explicit_unit'")
    if meta["value"] != "r":
        raise ValueError(f"{path}: metadata value must be 'r'")
    if not isinstance(meta["reference_checksum"], str) or meta["reference_checksum"] == "":
        raise ValueError(f"{path}: metadata reference_checksum must be a non-empty string")
    validate_plain_relative_filename(meta["reference_bim"], f"{path}: metadata reference_bim")
    if expected_format == NPZ_FORMAT:
        if meta.get("sparse_layout") != "csc":
            raise ValueError(f"{path}: metadata sparse_layout must be 'csc'")
        if meta.get("index_base") != 0:
            raise ValueError(f"{path}: metadata index_base must be 0")


def validate_chr_sex(chr_label, sex, where: str) -> None:
    if chr_label not in [str(i) for i in range(1, 23)] + ["X"]:
        raise ValueError(f"{where}: chr must be one of 1-22 or X")
    if chr_label == "X":
        validate_chrx_sex(sex, f"{where}: sex")
    elif sex is not None:
        raise ValueError(f"{where}: autosomal LD shards must have sex null")


def validate_chrx_sex(sex, where: str) -> str:
    if sex not in VALID_CHRX_SEX:
        raise ValueError(f"{where}: chrX sex must be one of female, male, combined")
    return sex


def validate_manifest_entry_agreement(entry: dict, meta: dict, path: Path) -> None:
    for key in ("chr", "sex", "num_snp", "nnz", "reference_checksum", "reference_bim"):
        if entry.get(key) != meta.get(key):
            raise ValueError(f"{path}: manifest/per-file metadata mismatch for {key}")
    suffix_runtime = runtime_from_suffix(path)
    if suffix_runtime == PY_RUNTIME_FORMAT and meta.get("format") != NPZ_FORMAT:
        raise ValueError(f"{path}: manifest file extension and metadata format disagree")


def default_reference_bim_filename(chr_label: str) -> str:
    return f"reference_chr{chr_label}.bim"


def validate_plain_relative_filename(value, where: str) -> str:
    if not isinstance(value, str) or value == "":
        raise ValueError(f"{where} must be a non-empty filename")
    path = Path(value)
    if path.is_absolute() or len(path.parts) != 1 or value in {".", ".."} or ".." in value:
        raise ValueError(f"{where} must be a plain relative filename")
    return value


validate_reference_bim_filename = validate_plain_relative_filename


def runtime_from_suffix(path: Path) -> str:
    if path.suffix == ".npz":
        return PY_RUNTIME_FORMAT
    if path.suffix == ".mat":
        raise ValueError(
            f"{path}: Python LD validation only supports .npz shards; "
            "use MATLAB/Octave to validate .mat LD distributions"
        )
    raise ValueError(f"{path}: unsupported LD shard extension {path.suffix!r}; expected .npz")


def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"LD file not found: {path}")


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_positive_int(value, where: str) -> None:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or int(value) <= 0:
        raise ValueError(f"{where} must be a positive integer")


def validate_nonnegative_int(value, where: str) -> None:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or int(value) < 0:
        raise ValueError(f"{where} must be a non-negative integer")
