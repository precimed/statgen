import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from statgen.ld import load_ld, validate_ld_distribution
from statgen.reference import load_reference
from tests.conftest import FIXTURES_DIR


SHARDED_BFILE = FIXTURES_DIR / "genotype/sharded/@"
SHARDED_REF = FIXTURES_DIR / "reference/sharded/@.bim"
BUILD_SCRIPT = Path(__file__).parents[1] / "script/statgen_build_ld.py"
MANIFEST_SCRIPT = Path(__file__).parents[1] / "script/statgen_create_ld_manifest.py"


def _metadata(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as z:
        return json.loads(bytes(z["metadata"]).decode("utf-8"))


def _write_fake_plink2(path: Path) -> None:
    path.write_text(
        r'''#!/usr/bin/env python3
import os
import sys
from pathlib import Path

args = sys.argv[1:]
if "--version" in args:
    print("PLINK v2.00 fake")
    raise SystemExit(0)

def arg(flag):
    if flag not in args:
        return None
    return args[args.index(flag) + 1]

log = os.environ.get("FAKE_PLINK_LOG")
if log:
    with open(log, "a") as f:
        f.write(" ".join(args) + "\n")

prefix = Path(arg("--bfile"))
out = Path(arg("--out"))
chr_label = arg("--chr")
keep_path = arg("--keep")
sample_count = 0
if keep_path:
    sample_count = sum(1 for line in open(keep_path) if line.strip())
else:
    sample_count = sum(1 for line in open(str(prefix) + ".fam") if line.strip())

rows = []
with open(str(prefix) + ".bim") as f:
    for raw in f:
        chrom, snp, cm, bp, a1, a2 = raw.split()
        if chr_label is None or chrom == chr_label:
            rows.append((chrom, snp, int(bp), a1, a2))

if "--freq" in args:
    with open(str(out) + ".afreq", "w") as f:
        f.write("#CHROM\tPOS\tID\tREF\tALT1\tALT1_FREQ\tOBS_CT\n")
        for i, (chrom, snp, bp, a1, a2) in enumerate(rows):
            freq = 0.0 if os.environ.get("FAKE_PLINK_MONOMORPHIC_FIRST") and i == 0 else 0.10 + 0.05 * i
            f.write(f"{chrom}\t{bp}\t{snp}\t{a2}\t{a1}\t{freq:.4f}\t{sample_count}\n")
    raise SystemExit(0)

if "--r-unphased" in args or "--r" in args:
    with open(str(out) + ".vcor", "w") as f:
        f.write("#CHROM_A\tPOS_A\tID_A\tREF_A\tALT1_A\tCHROM_B\tPOS_B\tID_B\tREF_B\tALT1_B\tUNPHASED_R\n")
        for i in range(max(0, len(rows) - 1)):
            a = rows[i]
            b = rows[i + 1]
            r = 0.25 if i % 2 == 0 else -0.40
            f.write(
                f"{a[0]}\t{a[2]}\t{a[1]}\t{a[4]}\t{a[3]}\t"
                f"{b[0]}\t{b[2]}\t{b[1]}\t{b[4]}\t{b[3]}\t{r}\n"
            )
    raise SystemExit(0)

print("unexpected fake plink invocation", args, file=sys.stderr)
raise SystemExit(2)
''',
        encoding="utf-8",
    )
    path.chmod(0o755)


def _load_build_module():
    from importlib.util import module_from_spec, spec_from_file_location

    spec = spec_from_file_location("statgen_build_ld", BUILD_SCRIPT)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _run_build(tmp_path: Path, args, *, expect_ok=True):
    fake = tmp_path / "fake_plink2.py"
    _write_fake_plink2(fake)
    log = tmp_path / "plink_commands.log"
    env = dict(**os_environ_without_pythonpath(), FAKE_PLINK_LOG=str(log))
    cmd = [sys.executable, str(BUILD_SCRIPT), "--plink2", str(fake), *args]
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if expect_ok:
        assert result.returncode == 0, result.stderr
    else:
        assert result.returncode != 0
    return result, log


def _create_manifest(out: Path):
    result = subprocess.run(
        [sys.executable, str(MANIFEST_SCRIPT), "--ld", str(out)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    return json.loads((out / "ld_manifest.json").read_text())


def _run_fake_build_for_shards(tmp_path: Path, out: Path, shards, extra_args=None):
    extra_args = [] if extra_args is None else list(extra_args)
    logs = []
    result = None
    for shard in shards:
        result, log = _run_build(
            tmp_path,
            [
                "--bfile",
                str(SHARDED_BFILE),
                "--out",
                str(out),
                "--shard",
                shard,
                *extra_args,
            ],
        )
        logs.append(log)
    return result, logs[-1]


def os_environ_without_pythonpath():
    import os

    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    return env


def _write_real_plink_mixed_vcf(path: Path) -> None:
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
        "1\t100\trs1\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1\n"
        "1\t200\trs2\tT\tC\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/0\n"
        "X\t3000000\tx1\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1\n"
        "X\t3000100\tx2\tT\tC\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/0\n",
        encoding="utf-8",
    )


def _write_real_plink_sex_psam(path: Path) -> None:
    path.write_text("#IID\tSEX\nS1\t1\nS2\t2\nS3\t1\nS4\t2\n", encoding="utf-8")


def _make_real_plink_bfile(plink2: str, tmp_path: Path, prefix_name: str, *, chr_filter=None) -> Path:
    vcf = tmp_path / f"{prefix_name}.vcf"
    psam = tmp_path / f"{prefix_name}.psam"
    _write_real_plink_mixed_vcf(vcf)
    _write_real_plink_sex_psam(psam)
    prefix = tmp_path / prefix_name
    cmd = [plink2, "--vcf", str(vcf), "--psam", str(psam), "--make-bed", "--out", str(prefix)]
    if chr_filter is not None:
        cmd.extend(["--chr", str(chr_filter)])
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return prefix


def test_statgen_build_ld_sharded_bfile_default_chrx_sex_split(tmp_path):
    out = tmp_path / "ld"
    _result, log = _run_fake_build_for_shards(tmp_path, out, ["1", "X"])
    manifest = _create_manifest(out)

    report = validate_ld_distribution(out, check_payload_structure=True)
    assert report["ok"] is True

    assert [(s["chr"], s["sex"]) for s in manifest["shards"]] == [
        ("1", None),
        ("X", "female"),
        ("X", "male"),
    ]
    assert [s["reference_bim"] for s in manifest["shards"]] == [
        "reference_chr1.bim",
        "reference_chrX.bim",
        "reference_chrX.bim",
    ]
    assert (out / "reference_chr1.bim").read_text() == (FIXTURES_DIR / "genotype/sharded/1.bim").read_text()
    assert (out / "reference_chrX.bim").read_text() == (FIXTURES_DIR / "genotype/sharded/X.bim").read_text()

    female_meta = _metadata(out / "ld_chrX_female.npz")
    male_meta = _metadata(out / "ld_chrX_male.npz")
    assert female_meta["reference_bim"] == "reference_chrX.bim"
    assert male_meta["reference_bim"] == "reference_chrX.bim"
    assert female_meta["num_sample"] == 2
    assert male_meta["num_sample"] == 2
    assert female_meta["plink_version"] == "PLINK v2.00 fake"
    assert "--r-unphased" in female_meta["build_command"]
    assert "--ld-window-kb 10000" in female_meta["build_command"]
    assert "--ld-window-r2 0.05" in female_meta["build_command"]

    commands = log.read_text()
    assert "--keep" in commands
    assert "--keep-allele-order" in commands

    reference = load_reference(SHARDED_REF)
    ld = load_ld(out)
    assert ld.reference.is_object_compatible(ld) is True
    assert [s.label for s in ld.shards] == ["1", "X"]
    assert ld.shard_groups[0][0].reference_checksum == reference.shards[0].checksum
    np.testing.assert_allclose(ld.a1freq("female")[:5], [0.10, 0.15, 0.20, 0.25, 0.30])
    assert ld.shard_groups[0][0].ld_r[0, 1] == pytest.approx(0.25)


def test_build_ld_distribution_python_api(tmp_path, monkeypatch):
    mod = _load_build_module()
    fake = tmp_path / "fake_plink2.py"
    _write_fake_plink2(fake)
    log = tmp_path / "api_plink_commands.log"
    monkeypatch.setenv("FAKE_PLINK_LOG", str(log))
    out = tmp_path / "ld_api"

    result = mod.build_ld_distribution(
        bfile=str(SHARDED_BFILE),
        out=out,
        shard="1",
        plink2=str(fake),
    )

    assert result is None
    assert not (out / "ld_manifest.json").exists()
    manifest = _create_manifest(out)
    assert [(record["chr"], record["sex"]) for record in manifest["shards"]] == [("1", None)]
    validate_ld_distribution(out, check_payload_structure=True)
    assert "--keep-allele-order" in log.read_text()
    assert (out / "_statgen_build_ld_work" / "chr1.afreq").is_file()


def test_create_ld_manifest_refuses_existing_manifest(tmp_path):
    out = tmp_path / "ld"
    _run_fake_build_for_shards(tmp_path, out, ["1"])
    _create_manifest(out)
    manifest_before = (out / "ld_manifest.json").read_text()

    result = subprocess.run(
        [sys.executable, str(MANIFEST_SCRIPT), "--ld", str(out)],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "Refusing to overwrite existing LD manifest" in result.stderr
    assert (out / "ld_manifest.json").read_text() == manifest_before


def test_build_ld_distribution_uses_explicit_scratch_directory(tmp_path, monkeypatch):
    mod = _load_build_module()
    fake = tmp_path / "fake_plink2.py"
    _write_fake_plink2(fake)
    log = tmp_path / "scratch_plink_commands.log"
    monkeypatch.setenv("FAKE_PLINK_LOG", str(log))
    out = tmp_path / "ld_scratch"
    scratch = tmp_path / "controlled_scratch"

    mod.build_ld_distribution(
        bfile=str(SHARDED_BFILE),
        out=out,
        shard="1",
        plink2=str(fake),
        scratch=scratch,
    )

    _create_manifest(out)
    validate_ld_distribution(out, check_payload_structure=True)
    assert (scratch / "chr1.afreq").is_file()
    assert not (out / "_statgen_build_ld_work").exists()
    assert str(scratch) in log.read_text()


def test_build_ld_distribution_rejects_negative_ld_r2_threshold(tmp_path):
    mod = _load_build_module()
    fake = tmp_path / "fake_plink2.py"
    _write_fake_plink2(fake)

    with pytest.raises(ValueError, match="ld_r2_threshold must be a finite non-negative value"):
        mod.build_ld_distribution(
            bfile=str(SHARDED_BFILE),
            out=tmp_path / "ld_negative_threshold",
            shard="1",
            plink2=str(fake),
            ld_r2_threshold=-0.01,
        )


def test_build_ld_distribution_rejects_monomorphic_snps_before_ld(tmp_path, monkeypatch):
    mod = _load_build_module()
    fake = tmp_path / "fake_plink2.py"
    _write_fake_plink2(fake)
    log = tmp_path / "monomorphic_reject.log"
    monkeypatch.setenv("FAKE_PLINK_LOG", str(log))
    monkeypatch.setenv("FAKE_PLINK_MONOMORPHIC_FIRST", "1")

    with pytest.raises(ValueError, match="contains 1 monomorphic SNPs"):
        mod.build_ld_distribution(
            bfile=str(SHARDED_BFILE),
            out=tmp_path / "ld_reject_monomorphic",
            shard="1",
            plink2=str(fake),
        )

    commands = log.read_text()
    assert "--freq" in commands
    assert "--r-unphased" not in commands
    assert not (tmp_path / "ld_reject_monomorphic" / "ld_chr1.monomorphic.tsv").exists()
    assert not (tmp_path / "ld_reject_monomorphic" / "reference_chr1.bim").exists()


def test_build_ld_distribution_allows_monomorphic_snps_with_sidecar(tmp_path, monkeypatch):
    mod = _load_build_module()
    fake = tmp_path / "fake_plink2.py"
    _write_fake_plink2(fake)
    log = tmp_path / "monomorphic_allow.log"
    monkeypatch.setenv("FAKE_PLINK_LOG", str(log))
    monkeypatch.setenv("FAKE_PLINK_MONOMORPHIC_FIRST", "1")
    out = tmp_path / "ld_allow_monomorphic"

    mod.build_ld_distribution(
        bfile=str(SHARDED_BFILE),
        out=out,
        shard="1",
        plink2=str(fake),
        allow_monomorphic_snps=True,
    )

    assert "--r-unphased" in log.read_text()
    sidecar = out / "ld_chr1.monomorphic.tsv"
    assert sidecar.read_text(encoding="utf-8").splitlines() == [
        "chr\tsnp\tbp\ta1\ta2",
        "1\trs1001\t100\tA\tG",
    ]
    assert _metadata(out / "ld_chr1.npz")["num_monomorphic_snps"] == 1


def test_statgen_build_ld_cli_passes_resource_controls(tmp_path):
    out = tmp_path / "ld"
    _result, log = _run_build(
        tmp_path,
        [
            "--bfile",
            str(SHARDED_BFILE),
            "--out",
            str(out),
            "--shard",
            "X",
            "--threads",
            "2",
            "--memory",
            "512",
        ],
    )

    commands = log.read_text()
    assert commands.count("--threads 2") == 4
    assert commands.count("--memory 512") == 4


def test_statgen_build_ld_cli_rejects_sample_missingness_filter(tmp_path):
    out = tmp_path / "ld"
    result, _log = _run_build(
        tmp_path,
        [
            "--bfile",
            str(SHARDED_BFILE),
            "--out",
            str(out),
            "--shard",
            "1",
            "--mind",
            "0.05",
        ],
        expect_ok=False,
    )

    assert "unrecognized arguments: --mind 0.05" in result.stderr


def test_build_ld_distribution_writes_chrx_split_keep_files_with_resource_controls(tmp_path, monkeypatch):
    mod = _load_build_module()
    fake = tmp_path / "fake_plink2.py"
    _write_fake_plink2(fake)
    log = tmp_path / "split_plink_commands.log"
    monkeypatch.setenv("FAKE_PLINK_LOG", str(log))
    out = tmp_path / "ld_split"
    scratch = tmp_path / "scratch_split"

    mod.build_ld_distribution(
        bfile=str(SHARDED_BFILE),
        out=out,
        shard="1",
        plink2=str(fake),
        scratch=scratch,
        threads=3,
        memory=1024,
    )
    mod.build_ld_distribution(
        bfile=str(SHARDED_BFILE),
        out=out,
        shard="X",
        plink2=str(fake),
        scratch=scratch,
        threads=3,
        memory=1024,
    )

    assert not (scratch / "chr1.keep").exists()
    assert (scratch / "chrX_male.keep").read_text() == "FAM1\tIND1\nFAM2\tIND3\n"
    assert (scratch / "chrX_female.keep").read_text() == "FAM1\tIND2\nFAM2\tIND4\n"
    chr1_meta = _metadata(out / "ld_chr1.npz")
    female_meta = _metadata(out / "ld_chrX_female.npz")
    male_meta = _metadata(out / "ld_chrX_male.npz")
    assert chr1_meta["num_sample"] == 4
    assert female_meta["num_sample"] == 2
    assert male_meta["num_sample"] == 2
    commands = log.read_text()
    assert commands.count("--threads 3") == 6
    assert commands.count("--memory 1024") == 6


def test_statgen_build_ld_end_to_end_with_real_plink2(tmp_path):
    plink2 = shutil.which("plink2")
    if plink2 is None:
        pytest.skip("plink2 not installed")

    prefix = _make_real_plink_bfile(plink2, tmp_path, "mixed")
    out = tmp_path / "ld_real"
    for shard in ["1", "X"]:
        build = subprocess.run(
            [
                sys.executable,
                str(BUILD_SCRIPT),
                "--plink2",
                plink2,
                "--bfile",
                str(prefix),
                "--out",
                str(out),
                "--shard",
                shard,
            ],
            capture_output=True,
            text=True,
        )
        assert build.returncode == 0, build.stderr
    _create_manifest(out)

    validate_ld_distribution(out, check_payload_structure=True)
    reference = load_reference(str(prefix) + ".bim")
    ld = load_ld(out)
    manifest = json.loads((out / "ld_manifest.json").read_text())
    assert [(s["chr"], s["sex"]) for s in manifest["shards"]] == [
        ("1", None),
        ("X", "female"),
        ("X", "male"),
    ]
    assert [s.label for s in ld.shards] == ["1", "X"]
    np.testing.assert_allclose(ld.a1freq("female"), [0.5, 0.375, 0.5, 0.25])
    assert ld.shards[0].ld_r.nnz >= 3
    np.testing.assert_allclose(ld.multiply_r2(np.ones(reference.num_snp)).shape, (4,))


def test_statgen_build_ld_real_plink2_combined_chrx(tmp_path):
    plink2 = shutil.which("plink2")
    if plink2 is None:
        pytest.skip("plink2 not installed")

    prefix = _make_real_plink_bfile(plink2, tmp_path, "mixed_combined")
    out = tmp_path / "ld_real_combined"
    for shard in ["1", "X"]:
        build = subprocess.run(
            [
                sys.executable,
                str(BUILD_SCRIPT),
                "--plink2",
                plink2,
                "--bfile",
                str(prefix),
                "--out",
                str(out),
                "--shard",
                shard,
                "--no-sex-split",
            ],
            capture_output=True,
            text=True,
        )
        assert build.returncode == 0, build.stderr
    _create_manifest(out)

    validate_ld_distribution(out, check_payload_structure=True)
    manifest = json.loads((out / "ld_manifest.json").read_text())
    assert [(s["chr"], s["sex"]) for s in manifest["shards"]] == [
        ("1", None),
        ("X", "combined"),
    ]
    x_meta = _metadata(out / "ld_chrX_combined.npz")
    assert "chrX_combined_rationale" not in x_meta
    ld = load_ld(out, default_chrX_sex="combined")
    assert ld.shards[1].sex == "combined"


def test_statgen_build_ld_real_plink2_sharded_input(tmp_path):
    plink2 = shutil.which("plink2")
    if plink2 is None:
        pytest.skip("plink2 not installed")

    _make_real_plink_bfile(plink2, tmp_path, "chr1", chr_filter="1")
    _make_real_plink_bfile(plink2, tmp_path, "chrX", chr_filter="X")
    out = tmp_path / "ld_real_sharded"
    for shard in ["1", "X"]:
        build = subprocess.run(
            [
                sys.executable,
                str(BUILD_SCRIPT),
                "--plink2",
                plink2,
                "--bfile",
                str(tmp_path / "chr@"),
                "--out",
                str(out),
                "--shard",
                shard,
            ],
            capture_output=True,
            text=True,
        )
        assert build.returncode == 0, build.stderr
    _create_manifest(out)

    validate_ld_distribution(out, check_payload_structure=True)
    manifest = json.loads((out / "ld_manifest.json").read_text())
    assert [(s["chr"], s["sex"]) for s in manifest["shards"]] == [
        ("1", None),
        ("X", "female"),
        ("X", "male"),
    ]
    ld = load_ld(out)
    assert [s.label for s in ld.shards] == ["1", "X"]


def test_statgen_build_ld_nonsharded_bfile_combined_chrx(tmp_path):
    prefix = tmp_path / "all"
    shutil.copyfile(FIXTURES_DIR / "reference/nonsharded/all.bim", str(prefix) + ".bim")
    shutil.copyfile(FIXTURES_DIR / "genotype/sharded/1.bed", str(prefix) + ".bed")
    shutil.copyfile(FIXTURES_DIR / "genotype/sharded/1.fam", str(prefix) + ".fam")
    out = tmp_path / "ld"

    for shard in ["1", "X"]:
        _run_build(
            tmp_path,
            [
                "--bfile",
                str(prefix),
                "--out",
                str(out),
                "--shard",
                shard,
                "--no-sex-split",
                "--ld-window-kb",
                "123",
                "--ld-r2-threshold",
                "0.2",
            ],
        )
    _create_manifest(out)

    validate_ld_distribution(out, check_payload_structure=True)
    manifest = json.loads((out / "ld_manifest.json").read_text())
    assert [(s["chr"], s["sex"]) for s in manifest["shards"]] == [
        ("1", None),
        ("X", "combined"),
    ]
    x_meta = _metadata(out / "ld_chrX_combined.npz")
    assert "chrX_combined_rationale" not in x_meta
    assert x_meta["ld_window_kb"] == 123
    assert x_meta["ld_r2_threshold"] == 0.2

    reference = load_reference(FIXTURES_DIR / "reference/nonsharded/all.bim")
    ld = load_ld(out, default_chrX_sex="combined")
    assert ld.shards[1].sex == "combined"
    np.testing.assert_allclose(ld.multiply_r2(np.ones(reference.num_snp))[:2], [1.0625, 1.2225])


def test_statgen_build_ld_rejects_missing_chrx_sex_without_no_sex_split(tmp_path):
    prefix = tmp_path / "X"
    shutil.copyfile(FIXTURES_DIR / "genotype/sharded/X.bim", str(prefix) + ".bim")
    shutil.copyfile(FIXTURES_DIR / "genotype/sharded/X.bed", str(prefix) + ".bed")
    (tmp_path / "X.fam").write_text(
        "FAM1\tIND1\t0\t0\t0\t-9\nFAM1\tIND2\t0\t0\t2\t-9\n",
        encoding="utf-8",
    )

    result, _log = _run_build(
        tmp_path,
        [
            "--bfile",
            str(prefix),
            "--out",
            str(tmp_path / "ld"),
            "--shard",
            "X",
        ],
        expect_ok=False,
    )
    assert "sex-specific builds require FAM sex codes 1 or 2" in (result.stderr + result.stdout)


def test_statgen_build_ld_rejects_empty_chrx_sex_group(tmp_path):
    prefix = tmp_path / "X"
    shutil.copyfile(FIXTURES_DIR / "genotype/sharded/X.bim", str(prefix) + ".bim")
    shutil.copyfile(FIXTURES_DIR / "genotype/sharded/X.bed", str(prefix) + ".bed")
    (tmp_path / "X.fam").write_text(
        "FAM1\tIND1\t0\t0\t1\t-9\nFAM1\tIND2\t0\t0\t1\t-9\n",
        encoding="utf-8",
    )

    result, _log = _run_build(
        tmp_path,
        [
            "--bfile",
            str(prefix),
            "--out",
            str(tmp_path / "ld"),
            "--shard",
            "X",
        ],
        expect_ok=False,
    )
    assert "chrX female build has no samples" in (result.stderr + result.stdout)


def test_parse_vcor_rejects_duplicate_unordered_pairs_with_path(tmp_path):
    from importlib.util import module_from_spec, spec_from_file_location

    spec = spec_from_file_location("statgen_build_ld", BUILD_SCRIPT)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)

    reference = load_reference(SHARDED_REF, shards=["1"])
    path = tmp_path / "dup.vcor"
    path.write_text(
        "#CHROM_A\tPOS_A\tREF_A\tALT1_A\tCHROM_B\tPOS_B\tREF_B\tALT1_B\tUNPHASED_R\n"
        "1\t100\tG\tA\t1\t200\tT\tC\t0.25\n"
        "1\t200\tT\tC\t1\t100\tG\tA\t0.25\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match=f"{path}: duplicate unordered LD pair"):
        mod._parse_vcor(path, reference.shards[0], chunk_rows=1)


def test_parse_vcor_accepts_chunked_table_reads(tmp_path):
    mod = _load_build_module()
    reference = load_reference(SHARDED_REF, shards=["1"])
    path = tmp_path / "chunked.vcor"
    path.write_text(
        "#CHROM_A\tPOS_A\tREF_A\tALT1_A\tCHROM_B\tPOS_B\tREF_B\tALT1_B\tUNPHASED_R\n"
        "1\t100\tG\tA\t1\t100\tG\tA\t1.0\n"
        "1\t100\tG\tA\t1\t200\tT\tC\t0.25\n"
        "1\t200\tT\tC\t1\t300\tC\tA\t-0.40\n",
        encoding="utf-8",
    )

    out = mod._parse_vcor(path, reference.shards[0], chunk_rows=1)

    assert list(out.columns) == ["idx1", "idx2", "r"]
    assert out["idx1"].to_numpy(dtype=np.int32).tolist() == [0, 1]
    assert out["idx2"].to_numpy(dtype=np.int32).tolist() == [1, 2]
    assert out["r"].to_numpy().dtype == np.float32
    np.testing.assert_allclose(out["r"], [0.25, -0.40])


def test_parse_afreq_rejects_missing_reference_variant_with_path(tmp_path):
    from importlib.util import module_from_spec, spec_from_file_location

    spec = spec_from_file_location("statgen_build_ld", BUILD_SCRIPT)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)

    reference = load_reference(SHARDED_REF, shards=["1"])
    path = tmp_path / "missing.afreq"
    path.write_text(
        "#CHROM\tPOS\tREF\tALT1\tALT1_FREQ\n"
        "1\t100\tG\tA\t0.10\n"
        "1\t200\tT\tC\t0.15\n"
        "1\t300\tC\tA\t0.20\n"
        "1\t400\tA\tG\t0.25\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match=f"{path}: missing frequency row for reference variant"):
        mod._parse_afreq(path, reference.shards[0])


def test_parse_vcor_rejects_absent_reference_variant_with_path(tmp_path):
    from importlib.util import module_from_spec, spec_from_file_location

    spec = spec_from_file_location("statgen_build_ld", BUILD_SCRIPT)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)

    reference = load_reference(SHARDED_REF, shards=["1"])
    path = tmp_path / "absent.vcor"
    path.write_text(
        "#CHROM_A\tPOS_A\tREF_A\tALT1_A\tCHROM_B\tPOS_B\tREF_B\tALT1_B\tUNPHASED_R\n"
        "1\t100\tG\tA\t1\t999\tT\tC\t0.25\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match=f"{path}: LD row references variant absent from reference shard"):
        mod._parse_vcor(path, reference.shards[0])


def test_plink2_missing_binary_and_nonzero_exit_errors(tmp_path):
    mod = _load_build_module()

    with pytest.raises(FileNotFoundError, match="PLINK2 executable not found"):
        mod._run([str(tmp_path / "missing_plink2")])

    failing = tmp_path / "failing_plink2.py"
    failing.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "print('plink failed intentionally', file=sys.stderr)\n"
        "raise SystemExit(1)\n",
        encoding="utf-8",
    )
    failing.chmod(0o755)

    with pytest.raises(RuntimeError, match="PLINK2 command failed:"):
        mod._run([str(failing), "--freq"])


def test_build_reference_rejects_reciprocal_snvs(tmp_path):
    mod = _load_build_module()
    prefix = tmp_path / "source"
    prefix.with_suffix(".bed").write_bytes(b"\x6c\x1b\x01\x00\x00")
    prefix.with_suffix(".fam").write_text("f\ti\t0\t0\t1\t-9\n")
    prefix.with_suffix(".bim").write_text(
        "1\ta\t0\t100\tA\tG\n1\tb\t0\t100\tG\tA\n")
    with pytest.raises(ValueError, match="reciprocal SNV"):
        mod._discover_bfile_shards(prefix, "1")


def test_ld_ordered_string_index_rejects_ambiguous_reference():
    from types import SimpleNamespace

    mod = _load_build_module()
    reference = SimpleNamespace(num_snp=2, chr=np.array(["1", "1"]),
                                bp=np.array([100, 100]), a1=np.array(["A", "A"]),
                                a2=np.array(["G", "G"]))
    with pytest.raises(ValueError, match="duplicate chr:bp:a1:a2"):
        mod._reference_key_index(reference)


def test_ld_output_matching_remains_exact_without_swaps(tmp_path):
    mod = _load_build_module()
    reference = load_reference(SHARDED_REF, shards=["1"]).shards[0]
    frequency = tmp_path / "reversed.afreq"
    frequency.write_text("#CHROM\tPOS\tREF\tALT1\tALT1_FREQ\n1\t100\tA\tG\t0.3\n")
    with pytest.raises(ValueError, match="frequency row references variant absent"):
        mod._parse_afreq(frequency, reference)
    correlations = tmp_path / "reversed.vcor"
    correlations.write_text(
        "#CHROM_A\tPOS_A\tREF_A\tALT1_A\tCHROM_B\tPOS_B\tREF_B\tALT1_B\tUNPHASED_R\n"
        "1\t100\tA\tG\t1\t200\tT\tC\t0.25\n")
    with pytest.raises(ValueError, match="LD row references variant absent"):
        mod._parse_vcor(correlations, reference)
