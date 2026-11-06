import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from statgen.genotype import load_genotype
from statgen.ld import load_ld, validate_ld_distribution
from tests.conftest import run_octave, skipif_no_octave


BUILD_SCRIPT = Path(__file__).parents[1] / "script/statgen_build_ld.py"
MANIFEST_SCRIPT = Path(__file__).parents[1] / "script/statgen_create_ld_manifest.py"


def _preflight_variant_rows(num_snp: int = 31):
    allele_pairs = [
        ("G", "A"),
        ("T", "C"),
        ("C", "G"),
        ("A", "T"),
        ("G", "C"),
        ("T", "A"),
    ]
    rows = []
    for i in range(num_snp):
        ref, alt = allele_pairs[i % len(allele_pairs)]
        rows.append(
            {
                "chrom": "1",
                "pos": 1000 + i * 100,
                "id": f"e2e_preflight_chr1_{i + 1:03d}_{alt}_{ref}",
                "ref": ref,
                "alt": alt,
            }
        )
    return rows


def _preflight_calls(num_snp: int = 31, num_sample: int = 29) -> np.ndarray:
    calls = np.fromfunction(
        lambda i, j: (i * 7 + j * 5 + (j % 4)) % 3,
        (num_snp, num_sample),
        dtype=int,
    ).astype(np.int8)

    sentinel = np.array(
        [0, 0, 0, 1, 1, 1, 2, 2, 2, 0, 1, 2, 0, 1, 2,
         0, 1, 2, 0, 1, 2, 0, 0, 1, 1, 2, 2, 0, 1],
        dtype=np.int8,
    )
    calls[0] = sentinel
    calls[1] = sentinel.copy()
    calls[1, [5, 17, 24]] = np.array([2, 0, 0], dtype=np.int8)
    calls[2] = 2 - sentinel
    calls[2, [8, 16, 25]] = np.array([0, 2, 1], dtype=np.int8)

    missing_positions = {
        1: [3, 21],
        2: [4, 20],
        5: [0, 11, 22],
        12: [2, 14],
        23: [6, 18, 27],
    }
    for snp_idx, sample_idx in missing_positions.items():
        calls[snp_idx, sample_idx] = -1

    for i, row in enumerate(calls):
        observed = row[row >= 0]
        if observed.size == 0 or np.unique(observed).size == 1:
            raise AssertionError(f"preflight generated monomorphic SNP row {i}")
    return calls


def _write_preflight_vcf(path: Path, variants, calls: np.ndarray) -> list[str]:
    samples = [f"S{i + 1:02d}" for i in range(calls.shape[1])]

    def gt(value: int) -> str:
        if value == -1:
            return "./."
        if value == 0:
            return "0/0"
        if value == 1:
            return "0/1"
        if value == 2:
            return "1/1"
        raise ValueError(value)

    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]
    for row, call_row in zip(variants, calls):
        lines.append(
            "\t".join(
                [
                    row["chrom"],
                    str(row["pos"]),
                    row["id"],
                    row["ref"],
                    row["alt"],
                    ".",
                    "PASS",
                    ".",
                    "GT",
                    *[gt(int(x)) for x in call_row],
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return samples


def _write_preflight_psam(path: Path, samples: list[str]) -> None:
    lines = ["#IID\tSEX"]
    for i, sample in enumerate(samples):
        if i < 14:
            sex = "1"
        elif i < 28:
            sex = "2"
        else:
            sex = "0"
        lines.append(f"{sample}\t{sex}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _run_checked(cmd: list[str]) -> subprocess.CompletedProcess:
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr + result.stdout
    return result


def _require_plink2() -> str:
    plink2 = shutil.which("plink2")
    if plink2 is None:
        pytest.skip("plink2 not installed")
    return plink2


def _create_manifest(out: Path) -> None:
    _run_checked([sys.executable, str(MANIFEST_SCRIPT), "--ld", str(out)])


def _read_bim(path: Path):
    rows = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        chrom, snp, cm, bp, a1, a2 = raw.split()
        rows.append({"chrom": chrom, "id": snp, "bp": int(bp), "a1": a1, "a2": a2})
    return rows


def _assert_bim_preserves_vcf_alleles(bim_path: Path, variants) -> None:
    bim = _read_bim(bim_path)
    assert len(bim) == len(variants)
    for got, expected in zip(bim, variants):
        assert got["chrom"] == expected["chrom"]
        assert got["id"] == expected["id"]
        assert got["bp"] == expected["pos"]
        assert got["a1"] == expected["alt"]
        assert got["a2"] == expected["ref"]


def _allele_frequency_from_fetch(genotypes: np.ndarray) -> np.ndarray:
    observed = np.isfinite(genotypes)
    denom = 2.0 * observed.sum(axis=0)
    return np.nansum(genotypes, axis=0) / denom


def _pearson_pair(genotypes: np.ndarray, i: int, j: int) -> float:
    x = genotypes[:, i]
    y = genotypes[:, j]
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 2:
        raise AssertionError("not enough observed calls for Pearson r")
    return float(np.corrcoef(x[keep], y[keep])[0, 1])


def _build_preflight_artifacts(tmp_path: Path) -> tuple[Path, Path]:
    plink2 = _require_plink2()
    variants = _preflight_variant_rows()
    calls = _preflight_calls()
    vcf = tmp_path / "preflight.vcf"
    psam = tmp_path / "preflight.psam"
    prefix = tmp_path / "preflight"
    samples = _write_preflight_vcf(vcf, variants, calls)
    _write_preflight_psam(psam, samples)

    _run_checked(
        [
            plink2,
            "--vcf",
            str(vcf),
            "--psam",
            str(psam),
            "--make-bed",
            "--keep-allele-order",
            "--out",
            str(prefix),
        ]
    )
    _assert_bim_preserves_vcf_alleles(prefix.with_suffix(".bim"), variants)

    ld_root = tmp_path / "ld"
    _run_checked(
        [
            sys.executable,
            str(BUILD_SCRIPT),
            "--plink2",
            plink2,
            "--bfile",
            str(prefix),
            "--out",
            str(ld_root),
            "--shard",
            "1",
        ]
    )
    _create_manifest(ld_root)
    assert validate_ld_distribution(ld_root, check_payload_structure=True)["ok"] is True
    return prefix, ld_root


def _python_preflight_outputs(prefix: Path, ld_root: Path) -> dict[str, np.ndarray]:
    ld = load_ld(ld_root)
    genotype = load_genotype(prefix, ld.reference)
    assert genotype.is_present.tolist() == [True] * ld.reference.num_snp

    fetched = genotype.fetch_genotypes(np.arange(genotype.num_snp))
    assert fetched.shape == (29, 31)
    ld_r = ld.shard_groups[0][0].ld_r
    vec = np.linspace(-1.5, 2.5, genotype.num_snp, dtype=np.float64)
    mat = np.column_stack([vec, vec[::-1]])
    return {
        "is_present": genotype.is_present.astype(np.float64),
        "ld_a1freq": ld.a1freq().astype(np.float64),
        "fetch_a1freq": _allele_frequency_from_fetch(fetched).astype(np.float64),
        "ld_r_selected": np.array([float(ld_r[0, 1]), float(ld_r[0, 2])], dtype=np.float64),
        "direct_r_selected": np.array(
            [_pearson_pair(fetched, 0, 1), _pearson_pair(fetched, 0, 2)],
            dtype=np.float64,
        ),
        "multiply_vec": ld.multiply_r2(vec).astype(np.float64),
        "explicit_vec": ((ld_r.toarray() ** 2) @ vec).astype(np.float64),
        "multiply_mat": ld.multiply_r2(mat).astype(np.float64).reshape(-1, order="F"),
        "explicit_mat": (((ld_r.toarray() ** 2) @ mat).astype(np.float64)).reshape(-1, order="F"),
    }


def _assert_preflight_outputs(outputs: dict[str, np.ndarray]) -> None:
    assert outputs["is_present"].tolist() == [1.0] * 31
    np.testing.assert_allclose(outputs["ld_a1freq"], outputs["fetch_a1freq"], atol=1e-5)
    for expected_r in outputs["direct_r_selected"]:
        assert abs(expected_r) > np.sqrt(0.05)
    np.testing.assert_allclose(outputs["ld_r_selected"], outputs["direct_r_selected"], atol=1e-3)
    np.testing.assert_allclose(
        outputs["multiply_vec"],
        outputs["explicit_vec"],
        atol=1e-5,
    )
    np.testing.assert_allclose(
        outputs["multiply_mat"],
        outputs["explicit_mat"],
        atol=1e-5,
    )


def _matlab_quote(path: Path) -> str:
    return str(path).replace("'", "''")


def _read_labeled_numeric_output(path: Path) -> dict[str, np.ndarray]:
    out = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        fields = raw.strip().split()
        if not fields:
            continue
        out[fields[0]] = np.array([float(x) for x in fields[1:]], dtype=np.float64)
    return out


def _full_sample_ids() -> list[str]:
    return [f"S{i + 1:02d}" for i in range(40)]


def _full_sample_sex(samples: list[str]) -> dict[str, str]:
    return {sample: ("1" if i < 20 else "2") for i, sample in enumerate(samples)}


def _full_x_subset(samples: list[str]) -> list[str]:
    male = samples[:20]
    female = samples[20:]
    return [
        male[18], female[19], male[0], female[1], male[5], female[6],
        male[12], female[13], male[2], female[3], male[8], female[9],
        male[14], female[15], male[4], female[5], male[10], female[11],
        male[16], female[17], male[6], female[7], male[19], female[18],
    ]


def _full_variant_rows():
    configs = {
        "1": {"n": 65, "shared": 32, "geno": 16, "ld": 17, "base": 1000},
        "2": {"n": 63, "shared": 31, "geno": 16, "ld": 16, "base": 2000},
        "X": {"n": 72, "shared": 37, "geno": 18, "ld": 17, "base": 3_000_000},
    }
    allele_pairs = [
        ("A", "C"),
        ("A", "G"),
        ("A", "T"),
        ("C", "A"),
        ("C", "G"),
        ("C", "T"),
        ("G", "A"),
        ("G", "C"),
        ("G", "T"),
        ("T", "A"),
        ("T", "C"),
        ("T", "G"),
    ]
    rows = []
    for chrom, cfg in configs.items():
        for i in range(cfg["n"]):
            if i < cfg["shared"]:
                category = "shared"
            elif i < cfg["shared"] + cfg["geno"]:
                category = "geno"
            else:
                category = "ld"
            ref, alt = allele_pairs[i % len(allele_pairs)]
            # The same-position cluster shares REF=A. Reciprocal SNVs at
            # one coordinate are invalid references under the swap contract.
            if i < 3:
                pos = cfg["base"]
            else:
                pos = cfg["base"] + 100 * (i - 3 + 1)
            rows.append(
                {
                    "chrom": chrom,
                    "pos": pos,
                    "id": f"e2e_full_chr{chrom}_{i + 1:03d}_{category}_{alt}_{ref}",
                    "ref": ref,
                    "alt": alt,
                    "category": category,
                    "ordinal": i,
                }
            )
    return rows


def _full_calls(variants, samples: list[str]) -> np.ndarray:
    calls = np.fromfunction(
        lambda i, j: (i * 5 + j * 7 + (j % 5)) % 3,
        (len(variants), len(samples)),
        dtype=int,
    ).astype(np.int8)

    base = np.array(([0, 0, 1, 1, 2, 2, 0, 1, 2, 0] * 4)[: len(samples)], dtype=np.int8)
    sample_index = {sample: i for i, sample in enumerate(samples)}
    by_id = {row["id"]: i for i, row in enumerate(variants)}

    for chrom in ("1", "2"):
        ids = [row["id"] for row in variants if row["chrom"] == chrom and row["category"] == "shared"]
        calls[by_id[ids[0]]] = base
        calls[by_id[ids[1]]] = base.copy()
        calls[by_id[ids[1]], [4, 17, 29]] = np.array([2, 0, 1], dtype=np.int8)
        calls[by_id[ids[2]]] = 2 - base
        calls[by_id[ids[2]], [9, 23, 35]] = np.array([1, 2, 0], dtype=np.int8)

    x_ids = [row["id"] for row in variants if row["chrom"] == "X" and row["category"] == "shared"]
    x0 = np.zeros(len(samples), dtype=np.int8)
    x1 = np.zeros(len(samples), dtype=np.int8)
    x2 = np.zeros(len(samples), dtype=np.int8)
    for j, sample in enumerate(samples):
        if j < 20:
            x0[j] = j % 2
            x1[j] = 1 - (j % 2)
            x2[j] = (j // 2) % 2
        else:
            val = [0, 1, 2][j % 3]
            x0[j] = val
            x1[j] = val
            x2[j] = 2 - val
    x1[[22, 31, 38]] = np.array([2, 0, 1], dtype=np.int8)
    calls[by_id[x_ids[0]]] = x0
    calls[by_id[x_ids[1]]] = x1
    calls[by_id[x_ids[2]]] = x2

    missing_specs = [
        (7, ["S03", "S18"]),
        (19, ["S11", "S27", "S36"]),
        (70, ["S02", "S21"]),
        (112, ["S16", "S34"]),
        (150, ["S05", "S24", "S39"]),
        (188, ["S09", "S29"]),
    ]
    for row_idx, sample_ids in missing_specs:
        for sample in sample_ids:
            calls[row_idx, sample_index[sample]] = -1

    for i, row in enumerate(calls):
        observed = row[row >= 0]
        if observed.size == 0 or np.unique(observed).size == 1:
            raise AssertionError(f"full e2e generated monomorphic SNP row {i}")
    return calls


def _write_psam(path: Path, samples: list[str], sex_by_sample: dict[str, str]) -> None:
    lines = ["#IID\tSEX"]
    lines.extend(f"{sample}\t{sex_by_sample[sample]}" for sample in samples)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_vcf(path: Path, variants, calls: np.ndarray, samples: list[str], all_samples: list[str]) -> None:
    sample_index = {sample: i for i, sample in enumerate(all_samples)}

    def gt(value: int, chrom: str, sample: str) -> str:
        if value == -1:
            return "."
        if chrom == "X" and sample_index[sample] < 20:
            return "1" if value > 0 else "0"
        if value == 0:
            return "0/0"
        if value == 1:
            return "0/1"
        if value == 2:
            return "1/1"
        raise ValueError(value)

    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]
    for row_idx, row in enumerate(variants):
        lines.append(
            "\t".join(
                [
                    row["chrom"],
                    str(row["pos"]),
                    row["id"],
                    row["ref"],
                    row["alt"],
                    ".",
                    "PASS",
                    ".",
                    "GT",
                    *[
                        gt(int(calls[row_idx, sample_index[sample]]), row["chrom"], sample)
                        for sample in samples
                    ],
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_id_file(path: Path, variants) -> None:
    path.write_text("".join(f"{row['id']}\n" for row in variants), encoding="utf-8")


def _plink_make_bed(plink2: str, vcf: Path, psam: Path, out_prefix: Path) -> None:
    _run_checked(
        [
            plink2,
            "--vcf",
            str(vcf),
            "--psam",
            str(psam),
            "--make-bed",
            "--keep-allele-order",
            "--out",
            str(out_prefix),
        ]
    )


def _plink_extract_bed(plink2: str, source_prefix: Path, extract: Path, out_prefix: Path, chrom=None) -> None:
    cmd = [
        plink2,
        "--bfile",
        str(source_prefix),
        "--extract",
        str(extract),
        "--make-bed",
        "--keep-allele-order",
        "--out",
        str(out_prefix),
    ]
    if chrom is not None:
        cmd.extend(["--chr", str(chrom)])
    _run_checked(cmd)


def _write_ploidy_for_bim(prefix: Path) -> None:
    rows = _read_bim(prefix.with_suffix(".bim"))
    values = [("1", "2") if row["chrom"] == "X" else ("2", "2") for row in rows]
    prefix.with_suffix(".ploidy").write_text(
        "".join(f"{male}\t{female}\n" for male, female in values),
        encoding="utf-8",
    )


def _build_ld_for_shards(plink2: str, bfile_prefix: Path, out_root: Path, shards=("1", "2", "X")) -> None:
    for shard in shards:
        _run_checked(
            [
                sys.executable,
                str(BUILD_SCRIPT),
                "--plink2",
                plink2,
                "--bfile",
                str(bfile_prefix),
                "--out",
                str(out_root),
                "--shard",
                shard,
            ]
        )
    _create_manifest(out_root)
    assert validate_ld_distribution(out_root, check_payload_structure=True)["ok"] is True


def _build_full_artifacts(tmp_path: Path) -> dict[str, Path]:
    plink2 = _require_plink2()
    variants = _full_variant_rows()
    samples = _full_sample_ids()
    sex_by_sample = _full_sample_sex(samples)
    x_samples = _full_x_subset(samples)
    calls = _full_calls(variants, samples)

    master_vcf = tmp_path / "master.vcf"
    master_psam = tmp_path / "master.psam"
    master_prefix = tmp_path / "master"
    _write_vcf(master_vcf, variants, calls, samples, samples)
    _write_psam(master_psam, samples, sex_by_sample)
    _plink_make_bed(plink2, master_vcf, master_psam, master_prefix)
    _assert_bim_preserves_vcf_alleles(master_prefix.with_suffix(".bim"), variants)

    genotype_variants = [row for row in variants if row["category"] in {"shared", "geno"}]
    ld_variants = [row for row in variants if row["category"] in {"shared", "ld"}]
    shared_variants = [row for row in variants if row["category"] == "shared"]
    genotype_ids = tmp_path / "genotype.extract"
    ld_ids = tmp_path / "ld.extract"
    _write_id_file(genotype_ids, genotype_variants)
    _write_id_file(ld_ids, ld_variants)

    genotype_prefix = tmp_path / "genotype_all"
    ld_source_prefix = tmp_path / "ld_source_all"
    _plink_extract_bed(plink2, master_prefix, genotype_ids, genotype_prefix)
    _plink_extract_bed(plink2, master_prefix, ld_ids, ld_source_prefix)
    _write_ploidy_for_bim(genotype_prefix)

    genotype_sharded = tmp_path / "genotype_sharded" / "chr@"
    ld_sharded = tmp_path / "ld_source_sharded" / "chr@"
    genotype_sharded.parent.mkdir()
    ld_sharded.parent.mkdir()
    for chrom in ("1", "2"):
        _plink_extract_bed(plink2, master_prefix, genotype_ids, tmp_path / "genotype_sharded" / f"chr{chrom}", chrom)
        _plink_extract_bed(plink2, master_prefix, ld_ids, tmp_path / "ld_source_sharded" / f"chr{chrom}", chrom)

    for stem, subset in [("genotype_sharded", genotype_variants), ("ld_source_sharded", ld_variants)]:
        x_variants = [row for row in subset if row["chrom"] == "X"]
        source_indices = [variants.index(row) for row in x_variants]
        x_vcf = tmp_path / f"{stem}_X.vcf"
        x_psam = tmp_path / f"{stem}_X.psam"
        x_prefix = tmp_path / stem / "chrX"
        _write_vcf(x_vcf, x_variants, calls[source_indices, :], x_samples, samples)
        _write_psam(x_psam, x_samples, sex_by_sample)
        _plink_make_bed(plink2, x_vcf, x_psam, x_prefix)
        _assert_bim_preserves_vcf_alleles(x_prefix.with_suffix(".bim"), x_variants)
        if stem == "genotype_sharded":
            _write_ploidy_for_bim(x_prefix)

    ld_full_root = tmp_path / "ld_full"
    ld_sharded_root = tmp_path / "ld_sharded"
    _build_ld_for_shards(plink2, ld_source_prefix, ld_full_root)
    _build_ld_for_shards(plink2, ld_sharded, ld_sharded_root)

    assert len(genotype_variants) == 150
    assert len(ld_variants) == 150
    assert len(shared_variants) == 100
    return {
        "genotype_prefix": genotype_prefix,
        "genotype_sharded": genotype_sharded,
        "ld_full_root": ld_full_root,
        "ld_sharded_root": ld_sharded_root,
    }


def _shard_bounds(reference):
    bounds = []
    start = 0
    for shard in reference.shards:
        stop = start + shard.num_snp
        bounds.append((shard.label, start, stop))
        start = stop
    return bounds


def _shard_for_index(bounds, index: int):
    for label, start, stop in bounds:
        if start <= index < stop:
            return label, index - start
    raise IndexError(index)


def _ploidy_adjusted_values(values: np.ndarray, ploidy: float) -> np.ndarray:
    if ploidy == 1:
        return np.minimum(values, 1.0)
    return values


def _expected_freq_from_genotype(genotype, chrX_sex: str) -> np.ndarray:
    present_indices = np.flatnonzero(genotype.is_present)
    fetched = genotype.fetch_genotypes(present_indices)
    out = np.full(genotype.num_snp, np.nan, dtype=np.float64)
    bounds = _shard_bounds(genotype)
    for col, global_idx in enumerate(present_indices):
        label, _local = _shard_for_index(bounds, int(global_idx))
        if label == "X":
            subject_mask = genotype.is_subject_present("X").copy()
            if chrX_sex == "female":
                subject_mask &= genotype.is_female
                ploidy = 2.0
            elif chrX_sex == "male":
                subject_mask &= genotype.is_male
                ploidy = 1.0
            else:
                raise ValueError(chrX_sex)
        else:
            subject_mask = np.ones(genotype.num_sample, dtype=bool)
            ploidy = 2.0
        values = fetched[subject_mask, col]
        observed = np.isfinite(values)
        if not observed.any():
            continue
        values = _ploidy_adjusted_values(values[observed], ploidy)
        out[global_idx] = values.sum() / (ploidy * observed.sum())
    return out


def _fetch_pair_for_ld(genotype, idx1: int, idx2: int, chrX_sex: str) -> tuple[np.ndarray, np.ndarray]:
    pair = genotype.fetch_genotypes([idx1, idx2])
    label1, _local1 = _shard_for_index(_shard_bounds(genotype), idx1)
    label2, _local2 = _shard_for_index(_shard_bounds(genotype), idx2)
    assert label1 == label2
    if label1 == "X":
        mask = genotype.is_subject_present("X").copy()
        if chrX_sex == "female":
            mask &= genotype.is_female
            ploidy = 2.0
        elif chrX_sex == "male":
            mask &= genotype.is_male
            ploidy = 1.0
        else:
            raise ValueError(chrX_sex)
    else:
        mask = np.ones(genotype.num_sample, dtype=bool)
        ploidy = 2.0
    x = pair[mask, 0]
    y = pair[mask, 1]
    keep = np.isfinite(x) & np.isfinite(y)
    return (
        _ploidy_adjusted_values(x[keep], ploidy),
        _ploidy_adjusted_values(y[keep], ploidy),
    )


def _direct_pair_r_for_ld(genotype, idx1: int, idx2: int, chrX_sex: str) -> float:
    x, y = _fetch_pair_for_ld(genotype, idx1, idx2, chrX_sex)
    if x.size < 2 or np.unique(x).size == 1 or np.unique(y).size == 1:
        raise AssertionError("selected pair is not variable enough for LD r check")
    return float(np.corrcoef(x, y)[0, 1])


def _find_global_index(reference, marker: str) -> int:
    matches = np.flatnonzero(reference.snp == marker)
    assert matches.size == 1, marker
    return int(matches[0])


def _ld_value(ld, idx1: int, idx2: int, chrX_sex: str) -> float:
    bounds = _shard_bounds(ld.reference)
    label1, local1 = _shard_for_index(bounds, idx1)
    label2, local2 = _shard_for_index(bounds, idx2)
    assert label1 == label2
    group = ld.shard_groups[[label for label, _start, _stop in bounds].index(label1)]
    if label1 == "X":
        shard = {s.sex: s for s in group}[chrX_sex]
    else:
        shard = group[0]
    return float(shard.ld_r[local1, local2])


def _full_selected_pair_indices(ld):
    pairs = []
    for chrom in ("1", "2", "X"):
        ids = [s for s in ld.reference.snp.tolist() if s.startswith(f"e2e_full_chr{chrom}_")]
        shared_ids = [s for s in ids if "_shared_" in s]
        pairs.append((_find_global_index(ld.reference, shared_ids[0]), _find_global_index(ld.reference, shared_ids[1])))
    return pairs


def _full_python_outputs(ld_root: Path, genotype_prefix: Path, *, sharded: bool) -> dict[str, np.ndarray]:
    ld = load_ld(ld_root)
    genotype = load_genotype(genotype_prefix, ld.reference)
    present = genotype.is_present
    assert int(present.sum()) == 100
    assert genotype.num_snp == 150

    freq_female = _expected_freq_from_genotype(genotype, "female")
    freq_male = _expected_freq_from_genotype(genotype, "male")
    present_indices = np.flatnonzero(present)
    pairs = _full_selected_pair_indices(ld)

    ld_r_values = []
    direct_r_values = []
    for pair in pairs[:2]:
        ld_r_values.append(_ld_value(ld, pair[0], pair[1], "female"))
        direct_r_values.append(_direct_pair_r_for_ld(genotype, pair[0], pair[1], "female"))
    x_pair = pairs[2]
    for sex in ("female", "male"):
        ld_r_values.append(_ld_value(ld, x_pair[0], x_pair[1], sex))
        direct_r_values.append(_direct_pair_r_for_ld(genotype, x_pair[0], x_pair[1], sex))

    vec = np.linspace(-2.0, 3.0, ld.num_snp, dtype=np.float64)
    multiply_female = ld.multiply_r2(vec, chrX_sex="female")
    multiply_male = ld.multiply_r2(vec, chrX_sex="male")

    outputs = {
        "present": present.astype(np.float64),
        "ld_a1freq_female_present": ld.a1freq("female")[present_indices].astype(np.float64),
        "calc_a1freq_female_present": freq_female[present_indices],
        "ld_a1freq_male_present": ld.a1freq("male")[present_indices].astype(np.float64),
        "calc_a1freq_male_present": freq_male[present_indices],
        "ld_r_selected": np.array(ld_r_values, dtype=np.float64),
        "direct_r_selected": np.array(direct_r_values, dtype=np.float64),
        "multiply_female_head": multiply_female[:10].astype(np.float64),
        "multiply_male_head": multiply_male[:10].astype(np.float64),
        "chrx_subject_present": genotype.is_subject_present("X").astype(np.float64),
    }
    if sharded:
        expected_mask = np.zeros(genotype.num_sample, dtype=np.float64)
        x_samples = set(_full_x_subset(_full_sample_ids()))
        for i, iid in enumerate(genotype.iid.tolist()):
            expected_mask[i] = 1.0 if iid in x_samples else 0.0
        np.testing.assert_array_equal(outputs["chrx_subject_present"], expected_mask)
    else:
        np.testing.assert_array_equal(outputs["chrx_subject_present"], np.ones(genotype.num_sample))
    return outputs


def _assert_full_outputs(outputs: dict[str, np.ndarray]) -> None:
    assert int(outputs["present"].sum()) == 100
    np.testing.assert_allclose(
        outputs["ld_a1freq_female_present"],
        outputs["calc_a1freq_female_present"],
        atol=1e-5,
    )
    np.testing.assert_allclose(
        outputs["ld_a1freq_male_present"],
        outputs["calc_a1freq_male_present"],
        atol=1e-5,
    )
    assert np.all(np.abs(outputs["direct_r_selected"]) > np.sqrt(0.05))
    np.testing.assert_allclose(outputs["ld_r_selected"], outputs["direct_r_selected"], atol=1e-3)


def test_preflight_chr1_a1freq_and_ld_r_match_fetched_genotypes(tmp_path):
    prefix, ld_root = _build_preflight_artifacts(tmp_path)
    _assert_preflight_outputs(_python_preflight_outputs(prefix, ld_root))


@pytest.mark.octave
@skipif_no_octave
def test_preflight_chr1_octave_matches_python_after_npz_to_mat_conversion(tmp_path):
    prefix, ld_root = _build_preflight_artifacts(tmp_path)
    expected = _python_preflight_outputs(prefix, ld_root)
    mat_root = tmp_path / "ld_matlab"
    out_path = tmp_path / "octave_preflight.tsv"

    script = (
        "warning('off', 'statgen:ld:v5mat'); "
        f"statgen.convert_ld_npz_to_mat('{_matlab_quote(ld_root)}', '{_matlab_quote(mat_root)}', '1', 'format', 'v5'); "
        f"report = statgen.create_ld_mat_manifest('{_matlab_quote(ld_root)}', '{_matlab_quote(mat_root)}', {{'1'}}); "
        f"validate_report = statgen.validate_ld_distribution('{_matlab_quote(mat_root)}', true); "
        f"ld = statgen.load_ld('{_matlab_quote(mat_root)}'); "
        f"g = statgen.load_genotype('{_matlab_quote(prefix)}', ld.reference); "
        "G = g.fetch_genotypes(1:g.num_snp); "
        "obs = isfinite(G); G0 = G; G0(~obs) = 0; "
        "freq = sum(G0, 1) ./ (2 * sum(obs, 1)); "
        "r01_keep = isfinite(G(:,1)) & isfinite(G(:,2)); "
        "r02_keep = isfinite(G(:,1)) & isfinite(G(:,3)); "
        "r01 = corrcoef(G(r01_keep,1), G(r01_keep,2)); "
        "r02 = corrcoef(G(r02_keep,1), G(r02_keep,3)); "
        "ld_r = ld.shard_groups{1}{1}.ld_r; "
        "v = linspace(-1.5, 2.5, g.num_snp)'; "
        "M = [v, flipud(v)]; "
        "mv = ld.multiply_r2(v); "
        "mm = ld.multiply_r2(M); "
        "explicit_v = full((ld_r .^ 2) * v); "
        "explicit_m = full((ld_r .^ 2) * M); "
        f"fid = fopen('{_matlab_quote(out_path)}', 'w'); "
        "fprintf(fid, 'manifest_count %.0f\\n', numel(report.shards)); "
        "fprintf(fid, 'validate_ok %.0f\\n', validate_report.ok); "
        "fprintf(fid, 'is_present '); fprintf(fid, '%.0f ', g.is_present); fprintf(fid, '\\n'); "
        "fprintf(fid, 'ld_a1freq '); fprintf(fid, '%.10g ', ld.a1freq()); fprintf(fid, '\\n'); "
        "fprintf(fid, 'fetch_a1freq '); fprintf(fid, '%.10g ', freq); fprintf(fid, '\\n'); "
        "fprintf(fid, 'ld_r_selected %.10g %.10g\\n', full(ld_r(1,2)), full(ld_r(1,3))); "
        "fprintf(fid, 'direct_r_selected %.10g %.10g\\n', r01(1,2), r02(1,2)); "
        "fprintf(fid, 'multiply_vec '); fprintf(fid, '%.10g ', mv); fprintf(fid, '\\n'); "
        "fprintf(fid, 'explicit_vec '); fprintf(fid, '%.10g ', explicit_v); fprintf(fid, '\\n'); "
        "fprintf(fid, 'multiply_mat '); fprintf(fid, '%.10g ', mm(:)); fprintf(fid, '\\n'); "
        "fprintf(fid, 'explicit_mat '); fprintf(fid, '%.10g ', explicit_m(:)); fprintf(fid, '\\n'); "
        "fclose(fid);"
    )
    result = run_octave(script, timeout=60)
    assert result.returncode == 0, result.stderr
    actual = _read_labeled_numeric_output(out_path)

    assert actual["manifest_count"].tolist() == [1.0]
    assert actual["validate_ok"].tolist() == [1.0]
    for key in [
        "is_present",
        "ld_a1freq",
        "fetch_a1freq",
        "ld_r_selected",
        "direct_r_selected",
        "multiply_vec",
        "explicit_vec",
        "multiply_mat",
        "explicit_mat",
    ]:
        np.testing.assert_allclose(actual[key], expected[key], atol=1e-5)

    _assert_preflight_outputs(actual)


def test_full_e2e_partial_overlap_a1freq_ld_r_and_chrx_subject_mapping_python(tmp_path):
    artifacts = _build_full_artifacts(tmp_path)

    nonsharded = _full_python_outputs(
        artifacts["ld_full_root"],
        artifacts["genotype_prefix"],
        sharded=False,
    )
    _assert_full_outputs(nonsharded)

    sharded = _full_python_outputs(
        artifacts["ld_sharded_root"],
        artifacts["genotype_sharded"],
        sharded=True,
    )
    _assert_full_outputs(sharded)

    assert sharded["chrx_subject_present"].sum() == 24
    assert nonsharded["chrx_subject_present"].sum() == 40


@pytest.mark.octave
@skipif_no_octave
def test_full_e2e_sharded_octave_matches_python_after_npz_to_mat_conversion(tmp_path):
    artifacts = _build_full_artifacts(tmp_path)
    expected = _full_python_outputs(
        artifacts["ld_sharded_root"],
        artifacts["genotype_sharded"],
        sharded=True,
    )
    ld_for_pairs = load_ld(artifacts["ld_sharded_root"])
    pairs = _full_selected_pair_indices(ld_for_pairs)
    bounds = _shard_bounds(ld_for_pairs.reference)
    pair_meta = []
    for idx1, idx2 in pairs:
        label, local1 = _shard_for_index(bounds, idx1)
        _label, local2 = _shard_for_index(bounds, idx2)
        group_idx = [b[0] for b in bounds].index(label) + 1
        pair_meta.append((idx1 + 1, idx2 + 1, group_idx, local1 + 1, local2 + 1))

    mat_root = tmp_path / "ld_sharded_matlab"
    out_path = tmp_path / "full_octave.tsv"

    p1 = pair_meta[0]
    p2 = pair_meta[1]
    px = pair_meta[2]
    script = (
        "warning('off', 'statgen:ld:v5mat'); "
        f"statgen.convert_ld_npz_to_mat('{_matlab_quote(artifacts['ld_sharded_root'])}', '{_matlab_quote(mat_root)}', '1', 'format', 'v5'); "
        f"statgen.convert_ld_npz_to_mat('{_matlab_quote(artifacts['ld_sharded_root'])}', '{_matlab_quote(mat_root)}', '2', 'format', 'v5'); "
        f"statgen.convert_ld_npz_to_mat('{_matlab_quote(artifacts['ld_sharded_root'])}', '{_matlab_quote(mat_root)}', 'X', 'format', 'v5'); "
        f"manifest = statgen.create_ld_mat_manifest('{_matlab_quote(artifacts['ld_sharded_root'])}', '{_matlab_quote(mat_root)}', {{'1', '2', 'X'}}); "
        f"validate_report = statgen.validate_ld_distribution('{_matlab_quote(mat_root)}', true); "
        f"ld = statgen.load_ld('{_matlab_quote(mat_root)}'); "
        f"g = statgen.load_genotype('{_matlab_quote(artifacts['genotype_sharded'])}', ld.reference); "
        "present_idx = find(g.is_present); "
        "G = g.fetch_genotypes(present_idx); "
        "freq_f = nan(g.num_snp, 1); freq_m = nan(g.num_snp, 1); "
        "n1 = ld.reference.shards{1}.num_snp; n2 = ld.reference.shards{2}.num_snp; x_start = n1 + n2 + 1; "
        "mask_x_f = g.is_subject_present('X') & g.is_female; "
        "mask_x_m = g.is_subject_present('X') & g.is_male; "
        "for c = 1:numel(present_idx); "
        "gi = present_idx(c); "
        "if gi >= x_start; "
        "vf = G(mask_x_f, c); of = isfinite(vf); freq_f(gi) = sum(vf(of)) / (2 * sum(of)); "
        "vm = G(mask_x_m, c); om = isfinite(vm); vm = min(vm(om), 1); freq_m(gi) = sum(vm) / sum(om); "
        "else; "
        "v = G(:, c); o = isfinite(v); f = sum(v(o)) / (2 * sum(o)); freq_f(gi) = f; freq_m(gi) = f; "
        "end; "
        "end; "
        f"P1 = g.fetch_genotypes([{p1[0]} {p1[1]}]); k = isfinite(P1(:,1)) & isfinite(P1(:,2)); C1 = corrcoef(P1(k,1), P1(k,2)); "
        f"P2 = g.fetch_genotypes([{p2[0]} {p2[1]}]); k = isfinite(P2(:,1)) & isfinite(P2(:,2)); C2 = corrcoef(P2(k,1), P2(k,2)); "
        f"PX = g.fetch_genotypes([{px[0]} {px[1]}]); "
        "xf = PX(mask_x_f,1); yf = PX(mask_x_f,2); kf = isfinite(xf) & isfinite(yf); CXF = corrcoef(xf(kf), yf(kf)); "
        "xm = min(PX(mask_x_m,1), 1); ym = min(PX(mask_x_m,2), 1); km = isfinite(xm) & isfinite(ym); CXM = corrcoef(xm(km), ym(km)); "
        "xgroup = ld.shard_groups{3}; "
        "if strcmp(xgroup{1}.sex, 'female'); xfemale = xgroup{1}; xmale = xgroup{2}; else; xfemale = xgroup{2}; xmale = xgroup{1}; end; "
        f"ld_r_selected = [full(ld.shard_groups{{{p1[2]}}}{{1}}.ld_r({p1[3]},{p1[4]})), "
        f"full(ld.shard_groups{{{p2[2]}}}{{1}}.ld_r({p2[3]},{p2[4]})), "
        f"full(xfemale.ld_r({px[3]},{px[4]})), full(xmale.ld_r({px[3]},{px[4]}))]; "
        "direct_r_selected = [C1(1,2), C2(1,2), CXF(1,2), CXM(1,2)]; "
        "v = linspace(-2.0, 3.0, ld.num_snp)'; "
        "mf = ld.multiply_r2(v, 'female'); mm = ld.multiply_r2(v, 'male'); "
        "af = ld.a1freq('female'); am = ld.a1freq('male'); "
        f"fid = fopen('{_matlab_quote(out_path)}', 'w'); "
        "fprintf(fid, 'manifest_count %.0f\\n', numel(manifest.shards)); "
        "fprintf(fid, 'validate_ok %.0f\\n', validate_report.ok); "
        "fprintf(fid, 'present '); fprintf(fid, '%.0f ', g.is_present); fprintf(fid, '\\n'); "
        "fprintf(fid, 'chrx_subject_present '); fprintf(fid, '%.0f ', g.is_subject_present('X')); fprintf(fid, '\\n'); "
        "fprintf(fid, 'ld_a1freq_female_present '); fprintf(fid, '%.10g ', af(present_idx)); fprintf(fid, '\\n'); "
        "fprintf(fid, 'calc_a1freq_female_present '); fprintf(fid, '%.10g ', freq_f(present_idx)); fprintf(fid, '\\n'); "
        "fprintf(fid, 'ld_a1freq_male_present '); fprintf(fid, '%.10g ', am(present_idx)); fprintf(fid, '\\n'); "
        "fprintf(fid, 'calc_a1freq_male_present '); fprintf(fid, '%.10g ', freq_m(present_idx)); fprintf(fid, '\\n'); "
        "fprintf(fid, 'ld_r_selected '); fprintf(fid, '%.10g ', ld_r_selected); fprintf(fid, '\\n'); "
        "fprintf(fid, 'direct_r_selected '); fprintf(fid, '%.10g ', direct_r_selected); fprintf(fid, '\\n'); "
        "fprintf(fid, 'multiply_female_head '); fprintf(fid, '%.10g ', mf(1:10)); fprintf(fid, '\\n'); "
        "fprintf(fid, 'multiply_male_head '); fprintf(fid, '%.10g ', mm(1:10)); fprintf(fid, '\\n'); "
        "fclose(fid);"
    )
    result = run_octave(script, timeout=120)
    assert result.returncode == 0, result.stderr
    actual = _read_labeled_numeric_output(out_path)

    assert actual["manifest_count"].tolist() == [4.0]
    assert actual["validate_ok"].tolist() == [1.0]
    for key in [
        "present",
        "chrx_subject_present",
        "ld_a1freq_female_present",
        "calc_a1freq_female_present",
        "ld_a1freq_male_present",
        "calc_a1freq_male_present",
        "ld_r_selected",
        "direct_r_selected",
        "multiply_female_head",
        "multiply_male_head",
    ]:
        np.testing.assert_allclose(actual[key], expected[key], atol=1e-5)
    _assert_full_outputs(actual)
