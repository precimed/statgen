#!/usr/bin/env python3
"""Shared helpers for current-schema Python cache fixture generators."""

from pathlib import Path
import sys

CACHE_DIR = Path(__file__).resolve().parent
FIXTURES_DIR = CACHE_DIR.parent
FIXTURES_REL = Path("tests/fixtures")
REPO_ROOT = FIXTURES_DIR.parent.parent
PYTHON_ROOT = REPO_ROOT / "python"
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

from statgen.annotations import load_annotation, load_annotations, save_annotations_cache
from statgen.genotype import load_genotype, save_genotype_cache
from statgen.reference import load_reference, save_reference_cache
from statgen.sumstats import load_sumstats, save_sumstats_cache


def _reference():
    return load_reference(FIXTURES_REL / "reference/sharded/@.bim")


def _annotations(reference):
    binary = load_annotations(
        [
            FIXTURES_REL / "annotations/anno1.bed",
            FIXTURES_REL / "annotations/anno2.bed",
        ],
        reference,
    )
    continuous = load_annotation(
        FIXTURES_REL / "annotations/continuous.annot",
        reference,
        has_header=True,
        value_columns=["score", "weight"],
        annotation_metadata_path=FIXTURES_REL / "annotations/continuous.meta",
    )
    return binary.union_annotations(continuous)


def generate_python_cache_fixture(object_name: str) -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    reference = (load_reference(FIXTURES_REL / "allele_swaps/reference.bim")
                 if object_name == "genotype" else _reference())
    if object_name == "reference":
        save_reference_cache(reference, CACHE_DIR / "reference_python_reference_cache_0_2.npz")
    elif object_name == "sumstats":
        panel = load_sumstats(FIXTURES_REL / "sumstats/traits_complete.tsv.gz", reference)
        save_sumstats_cache(panel, CACHE_DIR / "sumstats_python_sumstats_cache_0_1.npz")
    elif object_name == "annotations":
        save_annotations_cache(_annotations(reference), CACHE_DIR / "annotations_python_annotations_cache_0_2.npz")
    elif object_name == "genotype":
        panel = load_genotype(FIXTURES_REL / "allele_swaps/genotype/@", reference, allow_indel_swaps=True)
        save_genotype_cache(panel, CACHE_DIR / "genotype_python_genotype_cache_0_2.npz")
    else:
        raise ValueError(f"unknown cache fixture object: {object_name}")
