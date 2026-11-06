#!/usr/bin/env python3
"""Generate tests/fixtures/cache/genotype_python_genotype_cache_0_1.npz.

Run from the repository root:
    python tests/fixtures/cache/generate_genotype_python_genotype_cache_0_1.py
"""

from cache_fixture_helpers import generate_python_cache_fixture

from statgen.genotype import _CACHE_SCHEMA

if _CACHE_SCHEMA != "genotype_cache/0.1":
    raise SystemExit("Historical genotype /0.1 fixtures require revision 86e1656; see cache/README.md")

generate_python_cache_fixture("genotype")
