#!/usr/bin/env python3
"""
Generate all committed test fixtures under tests/fixtures/.

Run from anywhere:
    python tests/fixtures/generate.py

All binary arrays are written with explicit little-endian dtypes.
Reference checksums are MD5 over "chr:bp:a1:a2\n" lines in shard row order.
"""

import gzip
import hashlib
import json
import struct
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent


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


def write_plink_bed(path: Path, n_snp: int, n_sample: int) -> None:
    """Write a minimal PLINK .bed file (magic + SNP-major + all homref genotypes)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    bytes_per_snp = (n_sample + 3) // 4
    with open(path, "wb") as f:
        f.write(b"\x6c\x1b\x01")           # magic + SNP-major mode
        f.write(b"\x00" * (n_snp * bytes_per_snp))


def write_fam(path: Path, samples: list[tuple]) -> None:
    """Write a PLINK .fam file. samples = [(fid, iid, pid, mid, sex, pheno), ...]"""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for s in samples:
            f.write("\t".join(str(x) for x in s) + "\n")


def write_ld_shard(
    directory: Path,
    chr_label: str,
    sex: str | None,
    bim_rows: list[tuple],
    idx1: list[int],
    idx2: list[int],
    r: list[float],
    mafvec: list[float],
    num_sample: int = 100,
    extra_meta: dict | None = None,
) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    i1 = np.array(idx1, dtype="<i4")
    i2 = np.array(idx2, dtype="<i4")
    rv = np.array(r, dtype="<f4")
    maf = np.array(mafvec, dtype="<f4")
    i1.tofile(directory / "ld_idx1.i32")
    i2.tofile(directory / "ld_idx2.i32")
    rv.tofile(directory / "ld_r.f32")
    maf.tofile(directory / "mafvec.f32")
    meta = {
        "object_type": "ld_shard",
        "schema_version": "0.1",
        "chr": chr_label,
        "sex": sex,
        "num_snp": len(bim_rows),
        "num_ld": len(idx1),
        "index_base": 0,
        "byte_order": "little_endian",
        "triangle": "upper",
        "diagonal": "implicit_unit",
        "value": "r",
        "reference_checksum": bim_checksum(bim_rows),
        "build_tool": "generate.py",
        "build_command": "synthetic fixture",
        "ld_window_kb": 10000,
        "ld_r2_threshold": 0.05,
        "num_sample": num_sample,
    }
    if extra_meta:
        meta.update(extra_meta)
    with open(directory / "metadata.json", "w") as f:
        json.dump(meta, f, indent=2)
        f.write("\n")


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

# ---------------------------------------------------------------------------
# Fixture data
# ---------------------------------------------------------------------------

# LD for chr1 (upper-triangle pairs, idx1 < idx2)
CHR1_LD = dict(
    idx1=[0, 0, 1, 3],
    idx2=[1, 2, 3, 4],
    r=[0.9, 0.3, -0.5, 0.7],
    mafvec=[0.30, 0.40, 0.20, 0.35, 0.15],
)

# LD for chrX — same triplet structure, sex-specific mafvec and r
CHRX_LD_IDX1 = [0, 1]
CHRX_LD_IDX2 = [1, 2]
CHRX_LD = {
    "combined": dict(r=[0.60, -0.40], mafvec=[0.25, 0.30, 0.20]),
    "male":     dict(r=[0.55, -0.45], mafvec=[0.22, 0.28, 0.18]),
    "female":   dict(r=[0.65, -0.35], mafvec=[0.28, 0.32, 0.22]),
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

# Sumstats: chr, bp, a1, a2, z, n, p
# Row notes:
#   1:200  z=NA  (missing z)
#   1:300  p=0   (logp → Inf)
#   1:400  p=-0.1 (invalid p → logp NaN)
#   9:999  absent from reference → NaN row in aligned output
#   X:300 (rsX003) absent from file → NaN row in aligned output
SUMSTATS_ROWS = """\
chr\tbp\ta1\ta2\tz\tn\tp
1\t100\tA\tG\t2.5\t1000\t0.012
1\t200\tC\tT\tNA\t1000\t0.5
1\t300\tA\tC\t1.8\t1000\t0
1\t400\tG\tA\t-1.2\t1000\t-0.1
X\t100\tA\tG\t3.0\t500\t0.003
X\t200\tC\tT\t0.5\t500\t0.6
9\t999\tA\tG\t1.0\t1000\t0.3
"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    # --- reference ---
    write_bim(ROOT / "reference/sharded/1.bim", CHR1_BIM)
    write_bim(ROOT / "reference/sharded/X.bim", CHRX_BIM)
    write_bim(ROOT / "reference/nonsharded/all.bim", ALL_BIM)

    # --- annotations ---
    write_bed_annotation(ROOT / "annotations/anno1.bed", ANNO1_BED)
    write_bed_annotation(ROOT / "annotations/anno2.bed", ANNO2_BED)

    # --- sumstats ---
    sumstats_path = ROOT / "sumstats/traits.tsv.gz"
    sumstats_path.parent.mkdir(parents=True, exist_ok=True)
    # mtime=0 makes the gzip header deterministic (no embedded timestamp).
    with open(sumstats_path, "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
            gz.write(SUMSTATS_ROWS.encode())

    # --- LD chr1 ---
    write_ld_shard(
        ROOT / "ld/1",
        chr_label="1",
        sex=None,
        bim_rows=CHR1_BIM,
        **CHR1_LD,
    )

    # --- LD chrX shard group ---
    chrx_dir = ROOT / "ld/X"
    chrx_dir.mkdir(parents=True, exist_ok=True)
    sex_labels = list(CHRX_LD.keys())
    with open(chrx_dir / "metadata.json", "w") as f:
        json.dump(
            {"object_type": "ld_shard_group", "schema_version": "0.1",
             "chr": "X", "sex_shards": sex_labels},
            f, indent=2,
        )
        f.write("\n")
    for sex, data in CHRX_LD.items():
        write_ld_shard(
            chrx_dir / sex,
            chr_label="X",
            sex=sex,
            bim_rows=CHRX_BIM,
            idx1=CHRX_LD_IDX1,
            idx2=CHRX_LD_IDX2,
            **data,
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

    print("Fixtures written to", ROOT)
    print(f"  chr1 reference checksum : {bim_checksum(CHR1_BIM)}")
    print(f"  chrX reference checksum : {bim_checksum(CHRX_BIM)}")
    print(f"  all  reference checksum : {bim_checksum(ALL_BIM)}")


if __name__ == "__main__":
    main()
