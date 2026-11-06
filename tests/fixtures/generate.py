#!/usr/bin/env python3
"""
Generate all committed test fixtures under tests/fixtures/.

Run from anywhere:
    python tests/fixtures/generate.py

Binary fixture arrays use the dtypes documented by their object specs.
Reference checksums are MD5 over "chr:bp:a1:a2\n" lines in shard row order.
"""

import gzip
import hashlib
import io
import json
import shutil
import sys
import zipfile
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.io import savemat

ROOT = Path(__file__).parent
PYTHON_ROOT = ROOT.parents[1] / "python"
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from statgen._ld_schema import md5_file
from statgen._utils import allele_hash64
from statgen.reference import load_reference, save_reference_cache


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def bim_checksum(rows: list[tuple]) -> str:
    """MD5 over 'chr:bp:a1:a2\n' lines. rows = [(chr, snp, cm, bp, a1, a2), ...]"""
    text = "".join(f"{r[0]}:{r[3]}:{r[4]}:{r[5]}\n" for r in rows)
    return hashlib.md5(text.encode()).hexdigest()


def write_bim(path: Path, rows: list[tuple]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")


def write_bed_annotation(path: Path, intervals: list[tuple]) -> None:
    """Write a BED annotation file (chrom, start, end) — not PLINK .bed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for chrom, start, end in intervals:
            f.write(f"{chrom}\t{start}\t{end}\n")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_plink_bed(path: Path, n_snp: int, n_sample: int) -> None:
    """Write a minimal PLINK .bed file (magic + SNP-major + all homref genotypes)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    bytes_per_snp = (n_sample + 3) // 4
    with open(path, "wb") as f:
        f.write(b"\x6c\x1b\x01")           # magic + SNP-major mode
        f.write(b"\x00" * (n_snp * bytes_per_snp))


def write_ploidy(path: Path, rows: list[tuple]) -> None:
    """Write a .ploidy sidecar: tab-delimited, no header, (male_ploidy, female_ploidy) per row."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for male, female in rows:
            f.write(f"{male}\t{female}\n")


def write_fam(path: Path, samples: list[tuple]) -> None:
    """Write a PLINK .fam file. samples = [(fid, iid, pid, mid, sex, pheno), ...]"""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for s in samples:
            f.write("\t".join(str(x) for x in s) + "\n")


def ld_sparse_matrix(
    num_snp: int,
    idx1: list[int],
    idx2: list[int],
    r: list[float],
    dtype,
) -> sparse.csc_matrix:
    rows = np.array(idx1 + idx2 + list(range(num_snp)), dtype=np.int32)
    cols = np.array(idx2 + idx1 + list(range(num_snp)), dtype=np.int32)
    data = np.array(r + r + [1.0] * num_snp, dtype=dtype)
    return sparse.coo_matrix((data, (rows, cols)), shape=(num_snp, num_snp)).tocsc()


def ld_metadata(
    chr_label: str,
    sex: str | None,
    bim_rows: list[tuple],
    nnz: int,
    runtime_format: str,
    reference_bim: str | None = None,
    extra_meta: dict | None = None,
) -> dict:
    if reference_bim is None:
        reference_bim = f"reference_chr{chr_label}.bim"
    meta = {
        "object_type": "ld_shard",
        "schema_version": "1.0",
        "format": runtime_format,
        "chr": chr_label,
        "sex": sex,
        "num_snp": len(bim_rows),
        "nnz": nnz,
        "matrix": "symmetric",
        "diagonal": "explicit_unit",
        "value": "r",
        "num_monomorphic_snps": 0,
        "reference_checksum": bim_checksum(bim_rows),
        "reference_bim": reference_bim,
        "build_tool": "generate.py",
        "build_command": "synthetic fixture",
        "plink_version": None,
        "ld_window_kb": 10000,
        "ld_r2_threshold": 0.05,
        "num_sample": 100,
    }
    if runtime_format == "statgen_ld_npz_csc32":
        meta.update({"sparse_layout": "csc", "index_base": 0})
    if extra_meta:
        meta.update(extra_meta)
    return meta


def write_deterministic_npz(path: Path, arrays: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name, array in arrays.items():
            buf = io.BytesIO()
            np.save(buf, array, allow_pickle=False)
            info = zipfile.ZipInfo(f"{name}.npy", date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            zf.writestr(info, buf.getvalue())


def write_ld_npz_shard(
    path: Path,
    chr_label: str,
    sex: str | None,
    bim_rows: list[tuple],
    idx1: list[int],
    idx2: list[int],
    r: list[float],
    a1freq: list[float],
    extra_meta: dict | None = None,
) -> dict:
    mat = ld_sparse_matrix(len(bim_rows), idx1, idx2, r, np.float32)
    meta = ld_metadata(
        chr_label,
        sex,
        bim_rows,
        nnz=int(mat.nnz),
        runtime_format="statgen_ld_npz_csc32",
        extra_meta=extra_meta,
    )
    metadata_bytes = json.dumps(meta, sort_keys=True, separators=(",", ":")).encode()
    write_deterministic_npz(
        path,
        {
            "data": mat.data.astype(np.float32, copy=False),
            "indices": mat.indices.astype(np.int32, copy=False),
            "indptr": mat.indptr.astype(np.int32, copy=False),
            "shape": np.array(mat.shape, dtype=np.int64),
            "a1freq": np.array(a1freq, dtype=np.float32),
            "metadata": np.frombuffer(metadata_bytes, dtype=np.uint8),
        },
    )
    return meta


def matlab_metadata(meta: dict) -> dict:
    return {
        key: (np.array([], dtype=np.float64) if value is None else value)
        for key, value in meta.items()
    }


def write_ld_mat_shard(
    path: Path,
    chr_label: str,
    sex: str | None,
    bim_rows: list[tuple],
    idx1: list[int],
    idx2: list[int],
    r: list[float],
    a1freq: list[float],
    extra_meta: dict | None = None,
) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    mat = ld_sparse_matrix(len(bim_rows), idx1, idx2, r, np.float64)
    meta = ld_metadata(
        chr_label,
        sex,
        bim_rows,
        nnz=int(mat.nnz),
        runtime_format="statgen_ld_mat_sparse_double",
        extra_meta=extra_meta,
    )
    savemat(
        path,
        {
            "ld_r": mat,
            "a1freq": np.array(a1freq, dtype=np.float64).reshape(-1, 1),
            "metadata": matlab_metadata(meta),
        },
        appendmat=False,
        do_compression=False,
        long_field_names=True,
    )
    deterministic_header = (
        "MATLAB 5.0 MAT-file, Platform: statgen, Created by tests/fixtures/generate.py"
    )
    with open(path, "r+b") as f:
        f.write(deterministic_header.encode("ascii")[:116].ljust(116, b" "))
    return meta


def write_ld_manifest(
    directory: Path,
    runtime_format: str,
    shard_records: list[dict],
    reference_cache: str,
) -> None:
    manifest = {
        "object_type": "ld_panel_manifest",
        "schema_version": "1.0",
        "runtime_format": runtime_format,
        "reference_cache": reference_cache,
        "reference_cache_md5": md5_file(directory / reference_cache),
        "shards": shard_records,
    }
    with open(directory / "ld_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")


def write_ld_distribution(directory: Path, runtime_format: str, extension: str, shard_writer) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    reference_bims = {
        "1": "reference_chr1.bim",
        "X": "reference_chrX.bim",
    }
    write_bim(directory / reference_bims["1"], CHR1_BIM)
    write_bim(directory / reference_bims["X"], CHRX_BIM)
    if extension == "npz":
        reference_cache = "reference_cache.npz"
        reference = load_reference(str(directory / "reference_chr@.bim"))
        save_reference_cache(reference, directory / reference_cache)
    else:
        reference_cache = "reference_cache.mat"
        write_mat_reference_cache(directory / reference_cache, [("1", CHR1_BIM), ("X", CHRX_BIM)])
    shard_specs = [
        ("1", None, f"ld_chr1.{extension}", CHR1_BIM, CHR1_LD, None),
        ("X", "female", f"ld_chrX_female.{extension}", CHRX_BIM, CHRX_LD["female"], None),
        ("X", "male", f"ld_chrX_male.{extension}", CHRX_BIM, CHRX_LD["male"], None),
        (
            "X",
            "combined",
            f"ld_chrX_combined.{extension}",
            CHRX_BIM,
            CHRX_LD["combined"],
            {"chrX_combined_rationale": "synthetic fixture coverage"},
        ),
    ]
    records = []
    for chr_label, sex, file_name, bim_rows, data, extra_meta in shard_specs:
        path = directory / file_name
        idx1 = data.get("idx1", CHRX_LD_IDX1)
        idx2 = data.get("idx2", CHRX_LD_IDX2)
        meta = shard_writer(
            path,
            chr_label,
            sex,
            bim_rows,
            idx1,
            idx2,
            data["r"],
            data["a1freq"],
            extra_meta,
        )
        records.append({
            "chr": chr_label,
            "sex": sex,
            "file": file_name,
            "file_md5": md5_file(path),
            "num_snp": meta["num_snp"],
            "nnz": meta["nnz"],
            "reference_checksum": meta["reference_checksum"],
            "reference_bim": reference_bims[chr_label],
        })
    write_ld_manifest(directory, runtime_format, records, reference_cache)


def write_mat_reference_cache(path: Path, shards: list[tuple[str, list[tuple]]]) -> None:
    labels = [label for label, _rows in shards]
    checksums = [bim_checksum(rows) for _label, rows in shards]
    start0 = []
    stop0 = []
    pos = 0
    bp = []
    a1 = []
    a2 = []
    snp_text = []
    a1_text = []
    a2_text = []
    for _label, rows in shards:
        start0.append(pos)
        pos += len(rows)
        stop0.append(pos)
        bp.extend(int(r[3]) for r in rows)
        a1_values = [str(r[4]) for r in rows]
        a2_values = [str(r[5]) for r in rows]
        a1.extend(a1_values)
        a2.extend(a2_values)
        snp_text.append("\n".join(str(r[1]) for r in rows))
        a1_text.append("\n".join(a1_values))
        a2_text.append("\n".join(a2_values))

    metadata = {
        "schema": "reference_cache/0.1",
        "n_shards": len(shards),
        "shard_labels": np.array(labels, dtype=object).reshape(-1, 1),
        "shard_checksums": np.array(checksums, dtype=object).reshape(-1, 1),
        "shard_start0": np.array(start0, dtype=np.float64).reshape(-1, 1),
        "shard_stop0": np.array(stop0, dtype=np.float64).reshape(-1, 1),
    }
    savemat(
        path,
        {
            "metadata": metadata,
            "bp": np.array(bp, dtype=np.int64).reshape(-1, 1),
            "snp_text_by_shard": np.array(snp_text, dtype=object).reshape(-1, 1),
            "a1_text_by_shard": np.array(a1_text, dtype=object).reshape(-1, 1),
            "a2_text_by_shard": np.array(a2_text, dtype=object).reshape(-1, 1),
            "a1_hash64": allele_hash64(np.array(a1, dtype=object)).reshape(-1, 1),
            "a2_hash64": allele_hash64(np.array(a2, dtype=object)).reshape(-1, 1),
        },
        appendmat=False,
        do_compression=False,
        long_field_names=True,
    )
    deterministic_header = (
        "MATLAB 5.0 MAT-file, Platform: statgen, Created by tests/fixtures/generate.py"
    )
    with open(path, "r+b") as f:
        f.write(deterministic_header.encode("ascii")[:116].ljust(116, b" "))


# ---------------------------------------------------------------------------
# Variant definitions
# ---------------------------------------------------------------------------

# chr1: 5 SNPs — columns: chr, snp, cm, bp, a1, a2
CHR1_BIM = [
    ("1", "rs1001", 0, 100, "A", "G"),
    ("1", "rs1002", 0, 200, "C", "T"),
    ("1", "rs1003", 0, 300, "A", "C"),
    ("1", "rs1004", 0, 400, "G", "A"),
    ("1", "rs1005", 0, 500, "T", "C"),
]

# chrX: 3 SNPs
CHRX_BIM = [
    ("X", "rsX001", 0, 100, "A", "G"),
    ("X", "rsX002", 0, 200, "C", "T"),
    ("X", "rsX003", 0, 300, "A", "C"),
]

ALL_BIM = CHR1_BIM + CHRX_BIM

SAMPLES = [
    ("FAM1", "IND1", 0, 0, 1, -9),   # male
    ("FAM1", "IND2", 0, 0, 2, -9),   # female
    ("FAM2", "IND3", 0, 0, 1, -9),   # male
    ("FAM2", "IND4", 0, 0, 2, -9),   # female
]

# Ploidy: (male_ploidy, female_ploidy) per BIM row.
# chrX is hemizygous for males (1) and diploid for females (2).
CHRX_PLOIDY = [(1, 2)] * len(CHRX_BIM)
ALL_PLOIDY = [(2, 2)] * len(CHR1_BIM) + CHRX_PLOIDY

# ---------------------------------------------------------------------------
# Fixture data
# ---------------------------------------------------------------------------

# LD for chr1 (upper-triangle pairs, idx1 < idx2)
CHR1_LD = dict(
    idx1=[0, 0, 1, 3],
    idx2=[1, 2, 3, 4],
    r=[0.9, 0.3, -0.5, 0.7],
    a1freq=[0.30, 0.40, 0.20, 0.35, 0.15],
)

# LD for chrX — same sparse structure, sex-specific a1freq and r
CHRX_LD_IDX1 = [0, 1]
CHRX_LD_IDX2 = [1, 2]
CHRX_LD = {
    "female":   dict(r=[0.65, -0.35], a1freq=[0.28, 0.32, 0.22]),
    "male":     dict(r=[0.55, -0.45], a1freq=[0.22, 0.28, 0.18]),
    "combined": dict(r=[0.60, -0.40], a1freq=[0.25, 0.30, 0.20]),
}

# anno1.bed — chr1; overlapping [99,150)+[120,200) and adjacent [299,350)+[350,400)
#   SNP bp→0-based: 100→99, 200→199, 300→299, 400→399, 500→499
#   [99,150) ∪ [120,200) → after merge: [99,200)  covers rs1001(99), rs1002(199)
#   [299,350) + [350,400) → adjacent, merged: [299,400)  covers rs1003(299), rs1004(399)
#   [499,500) → boundary: covers rs1005(499) exactly
ANNO1_BED = [
    ("1", 99,  150),   # overlapping pair ↓
    ("1", 120, 200),
    ("1", 299, 350),   # adjacent pair ↓
    ("1", 350, 400),
    ("1", 499, 500),   # boundary: exactly covers rs1005
]

# anno2.bed — chrX; overlapping [195,205)+[198,302)
#   rsX001(99): covered by [99,100)
#   rsX002(199): covered by merged [195,302)
#   rsX003(299): covered by merged [195,302)
ANNO2_BED = [
    ("X",  99, 100),   # covers rsX001 only
    ("X", 195, 205),   # overlapping pair ↓
    ("X", 198, 302),
]

CONTINUOUS_ANNOTATION_ROWS = """\
chrom\tstart0\tend0\tscore\tweight
1\t99\t200\t0.5\t1.25
1\t299\t400\t-1.25\t0.75
X\t99\t301\t2.0\t3.5
"""

CONTINUOUS_ANNOTATION_META = """\
chrom metadata
start metadata
end metadata
score metadata
weight metadata
"""

GROUPED_ANNOTATION_ROWS = """\
chrom\tstart0\tend0\tgroup
1\t99\t200\tcoding
1\t150\t350\tregulatory
1\t299\t400\tcoding
X\t99\t301\tregulatory
"""

# Sumstats: chr, bp, a1, a2, z, n, p
# Row notes:
#   1:200  z=NA  (missing z)
#   1:300  p=0   (logp → Inf)
#   9:999  absent from reference → NaN row in aligned output
#   X:300 (rsX003) absent from file → NaN row in aligned output
SUMSTATS_ROWS = """\
chr\tbp\ta1\ta2\tz\tn\tp
1\t100\tA\tG\t2.5\t1000\t0.012
1\t200\tC\tT\tNA\t1000\t0.5
1\t300\tA\tC\t1.8\t1000\t0
1\t400\tG\tA\t-1.2\t1000\t0.2
X\t100\tA\tG\t3.0\t500\t0.003
X\t200\tC\tT\t0.5\t500\t0.6
9\t999\tA\tG\t1.0\t1000\t0.3
"""

# Complete, warning-free toy sumstats for package examples and vignettes.
SUMSTATS_EXAMPLE_ROWS = """\
chr\tbp\ta1\ta2\tz\tn\tp
1\t100\tA\tG\t2.5\t1000\t0.012
1\t200\tC\tT\t1.0\t1000\t0.5
1\t300\tA\tC\t1.8\t1000\t0.05
1\t400\tG\tA\t-1.2\t1000\t0.2
1\t500\tT\tC\t0.7\t1000\t0.4
X\t100\tA\tG\t3.0\t500\t0.003
X\t200\tC\tT\t0.5\t500\t0.6
X\t300\tG\tA\t-0.8\t500\t0.7
"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def write_allele_swap_fixtures() -> None:
    """Shared valid inputs and an oracle built from explicit construction choices.

    Mapping indices are zero-based physical rows; -1 denotes absence. JSON null
    denotes missing numeric data. Expected genotype matrices have sample rows
    and reference-variant columns, including null columns for absent variants
    (fetching an absent variant must itself fail). BIM/FAM are UTF-8 TSV; BED
    is standard SNP-major PLINK 1; expected arrays use JSON numeric/bool values.
    """
    root = ROOT / "allele_swaps"
    ref1 = [
        ("1", "exact", 0, 100, "A", "G"),
        ("1", "swap", 0, 200, "A", "G"),
        ("1", "amb_at", 0, 300, "A", "T"),
        ("1", "amb_cg", 0, 400, "C", "G"),
        ("1", "indel", 0, 500, "A", "AG"),
        ("1", "recip_indel_1", 0, 600, "A", "AT"),
        ("1", "recip_indel_2", 0, 600, "AT", "A"),
        ("1", "long_exact", 0, 700, "AC", "GT"),
        ("1", "long_no_swap", 0, 800, "AC", "GT"),
        ("1", "no_reuse_1", 0, 900, "A", "AC"),
        ("1", "no_reuse_2", 0, 900, "AC", "A"),
        ("1", "missing_exact", 0, 1000, "A", "G"),
        ("1", "absent", 0, 1100, "C", "T"),
    ]
    refx = [
        ("X", "x_swap", 0, 100, "A", "G"),
        ("X", "x_exact", 0, 200, "C", "T"),
        ("X", "x_indel", 0, 300, "A", "AT"),
    ]
    source1 = [
        ("1", "s0", 0, 100, "A", "G"),
        ("1", "s1", 0, 100, "G", "A"),
        ("1", "s2", 0, 200, "G", "A"),
        ("1", "s3", 0, 300, "T", "A"),
        ("1", "s4", 0, 400, "G", "C"),
        ("1", "s5", 0, 500, "AG", "A"),
        ("1", "s6", 0, 600, "A", "AT"),
        ("1", "s7", 0, 600, "AT", "A"),
        ("1", "s8", 0, 700, "AC", "GT"),
        ("1", "s9", 0, 700, "GT", "AC"),
        ("1", "s10", 0, 800, "GT", "AC"),
        ("1", "s11", 0, 900, "A", "AC"),
        ("1", "s12", 0, 1000, "A", "G"),
        ("1", "s13", 0, 1000, "G", "A"),
    ]
    sourcex = [
        ("X", "xs0", 0, 100, "G", "A"),
        ("X", "xs1", 0, 200, "C", "T"),
        ("X", "xs2", 0, 300, "AT", "A"),
    ]
    # Source SNP rows, subject columns. The exact row at 1000 is all missing;
    # its reciprocal row has complete data and must still lose.
    calls1 = [[0, 1, 2, -1], [2, 2, 2, 2], [0, 2, 1, -1],
              [2, 0, -1, 1], [1, 2, 0, -1], [0, 1, 2, -1],
              [2, 1, 0, -1], [0, 2, 1, -1], [1, 0, 2, -1],
              [2, 2, 2, 2], [0, 1, 2, -1], [2, 0, 1, -1],
              [-1, -1, -1, -1], [2, 2, 2, 2]]
    # X subjects are a reordered subset: sample 2 (male), then sample 1 (female).
    callsx = [[0, 2], [1, -1], [2, 0]]
    ploidyx = [(1, 2), (2, 2), (0, 1)]
    for label, rows, samples, calls in [
        ("1", source1, SAMPLES, calls1),
        ("X", sourcex, [SAMPLES[2], SAMPLES[1]], callsx),
    ]:
        base = root / f"genotype/{label}"
        write_bim(base.with_suffix(".bim"), rows)
        write_fam(base.with_suffix(".fam"), samples)
        # Explicit calls -> PLINK packed bytes, independent of library decoder.
        bits = {2: 0b00, -1: 0b01, 1: 0b10, 0: 0b11}
        packed = bytearray(b"\x6c\x1b\x01")
        for row in calls:
            for start in range(0, len(row), 4):
                packed.append(sum(bits[value] << (2 * j)
                                  for j, value in enumerate(row[start:start + 4])))
        base.with_suffix(".bed").write_bytes(packed)
    write_ploidy(root / "genotype/X.ploidy", ploidyx)
    write_bim(root / "reference.bim", ref1 + refx)

    sources = source1 + sourcex
    # p/n/se/info remain unchanged; missing optional values deliberately occur
    # in both an exact row and a swapped row.
    stats = [dict(p=0.01 * (j + 1), z=float(j + 1), n=100 + j,
                  beta=(j + 1) / 10, se=0.1, eaf=0.2, info=0.9)
             for j in range(len(sources))]
    for field in ("z", "n", "beta", "se", "eaf", "info"):
        stats[12][field] = None
    stats[3]["eaf"] = None
    columns = ["p", "z", "n", "beta", "se", "eaf", "info"]
    # Reverse TSV source order to exercise alignment independently of BIM order.
    tsv = "chr\tbp\ta1\ta2\t" + "\t".join(columns) + "\n"
    for row, values in reversed(list(zip(sources, stats))):
        tsv += "\t".join([row[0], str(row[3]), row[4], row[5]] +
                         ["NA" if values[c] is None else str(values[c]) for c in columns]) + "\n"
    with open(root / "sumstats.tsv.gz", "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
            gz.write(tsv.encode("utf-8"))
    expected = {"index_base": 0, "missing_source_row": -1,
                "reference_shard_sizes": {"1": 13, "X": 3}, "options": {}}
    for enabled in (False, True):
        # These selections and orientations are the independent oracle. No
        # matcher or loader is called to produce expected aligned outputs.
        mapping = [0, 2, 3, 4, 5 if enabled else -1, 6, 7, 8, -1,
                   11, -1, 12, -1, 14, 15, 16 if enabled else -1]
        swapped = [False, True, True, True, enabled, False, False, False,
                   False, False, False, False, False, True, False, enabled]
        aligned = {c: [] for c in columns}
        raw_columns, scaled_columns = [], []
        for i, (j, flip) in enumerate(zip(mapping, swapped)):
            for c in columns:
                value = None if j < 0 else stats[j][c]
                if value is not None and flip:
                    if c in ("z", "beta"):
                        value = -value
                    elif c == "eaf":
                        value = 1 - value
                aligned[c].append(value)
            if j < 0:
                raw = [None] * 4
            elif j < 14:
                raw = [None if v == -1 else v for v in calls1[j]]
            else:
                x = callsx[j - 14]
                raw = [None, None if x[1] == -1 else x[1],
                       None if x[0] == -1 else x[0], None]
            if flip:
                raw = [None if v is None else 2 - v for v in raw]
            male, female = (2, 2) if i < 13 else ploidyx[i - 13]
            scaled = [None if v is None else v * p / 2
                      for v, p in zip(raw, [male, female, male, female])]
            raw_columns.append(raw)
            scaled_columns.append(scaled)
        expected["options"][str(enabled).lower()] = {
            "source_row0": mapping[:13] + [j - 14 if j >= 0 else -1 for j in mapping[13:]],
            "sumstats_source_row0": [len(sources) - 1 - j if j >= 0 else -1 for j in mapping],
            "is_allele_swapped": swapped, "aligned_sumstats": aligned,
            "raw_counts": list(map(list, zip(*raw_columns))),
            "scaled_counts": list(map(list, zip(*scaled_columns))),
            "reporting": {"1": {"swaps": 4 if enabled else 3, "strand_ambiguous": 2,
                                  "disabled_indels": 0 if enabled else 1},
                          "X": {"swaps": 2 if enabled else 1, "strand_ambiguous": 0,
                                  "disabled_indels": 0 if enabled else 1}},
        }
    write_text(root / "expected.json", json.dumps(expected, indent=2, allow_nan=False) + "\n")


def main() -> None:
    write_allele_swap_fixtures()
    # --- reference ---
    write_bim(ROOT / "reference/sharded/1.bim", CHR1_BIM)
    write_bim(ROOT / "reference/sharded/X.bim", CHRX_BIM)
    write_bim(ROOT / "reference/nonsharded/all.bim", ALL_BIM)

    # --- annotations ---
    write_bed_annotation(ROOT / "annotations/anno1.bed", ANNO1_BED)
    write_bed_annotation(ROOT / "annotations/anno2.bed", ANNO2_BED)
    write_text(ROOT / "annotations/continuous.annot", CONTINUOUS_ANNOTATION_ROWS)
    write_text(ROOT / "annotations/continuous.meta", CONTINUOUS_ANNOTATION_META)
    write_text(ROOT / "annotations/grouped.annot", GROUPED_ANNOTATION_ROWS)

    # --- sumstats ---
    for name, rows in [
        ("traits.tsv.gz", SUMSTATS_ROWS),
        ("traits_complete.tsv.gz", SUMSTATS_EXAMPLE_ROWS),
    ]:
        sumstats_path = ROOT / f"sumstats/{name}"
        sumstats_path.parent.mkdir(parents=True, exist_ok=True)
        # mtime=0 makes the gzip header deterministic (no embedded timestamp).
        with open(sumstats_path, "wb") as raw:
            with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
                gz.write(rows.encode())

    # --- LD runtime distributions ---
    shutil.rmtree(ROOT / "ld", ignore_errors=True)
    write_ld_distribution(
        ROOT / "ld/python",
        runtime_format="python_npz_csc32",
        extension="npz",
        shard_writer=write_ld_npz_shard,
    )
    write_ld_distribution(
        ROOT / "ld/matlab",
        runtime_format="matlab_mat_sparse_double",
        extension="mat",
        shard_writer=write_ld_mat_shard,
    )

    # --- genotype (PLINK bfile metadata) ---
    for label, bim_rows in [("1", CHR1_BIM), ("X", CHRX_BIM)]:
        base = ROOT / f"genotype/sharded/{label}"
        write_bim(Path(str(base) + ".bim"), bim_rows)
        write_fam(Path(str(base) + ".fam"), SAMPLES)
        write_plink_bed(
            Path(str(base) + ".bed"),
            n_snp=len(bim_rows),
            n_sample=len(SAMPLES),
        )
    write_ploidy(ROOT / "genotype/sharded/X.ploidy", CHRX_PLOIDY)

    # --- genotype non-sharded bfile ---
    ns_base = ROOT / "genotype/nonsharded/all"
    write_bim(Path(str(ns_base) + ".bim"), ALL_BIM)
    write_fam(Path(str(ns_base) + ".fam"), SAMPLES)
    write_plink_bed(Path(str(ns_base) + ".bed"), n_snp=len(ALL_BIM), n_sample=len(SAMPLES))
    write_ploidy(Path(str(ns_base) + ".ploidy"), ALL_PLOIDY)

    print("Fixtures written to", ROOT)
    print(f"  chr1 reference checksum : {bim_checksum(CHR1_BIM)}")
    print(f"  chrX reference checksum : {bim_checksum(CHRX_BIM)}")
    print(f"  all  reference checksum : {bim_checksum(ALL_BIM)}")


if __name__ == "__main__":
    main()
