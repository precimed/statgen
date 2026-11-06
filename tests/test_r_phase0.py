"""R phase 0 acceptance tests: package skeleton and optional R smoke tests."""

import json
import os
import subprocess

import pytest

from tests.conftest import (
    FIXTURES_DIR,
    R_EXTDATA_COPIES,
    R_PACKAGE_DIR,
    REPO_ROOT,
    run_rscript,
    skipif_no_rscript,
)


def _r_package_version() -> str:
    for line in (R_PACKAGE_DIR / "DESCRIPTION").read_text().splitlines():
        if line.startswith("Version:"):
            return line.split(":", 1)[1].strip()
    raise AssertionError("R package DESCRIPTION is missing Version")


def test_r_package_skeleton_present():
    required = [
        "DESCRIPTION",
        "NAMESPACE",
        "R/version.R",
        "R/verbosity.R",
        "man/statgen-package.Rd",
        "man/version.Rd",
        "man/verbosity.Rd",
        "tests/testthat.R",
        "tests/testthat/test-ld.R",
        "tests/testthat/test-smoke.R",
        "inst/extdata/README.md",
    ]
    for rel in required:
        assert (R_PACKAGE_DIR / rel).is_file(), f"missing R package file: {rel}"


@pytest.mark.r
@skipif_no_rscript
def test_r_extdata_reference_fixture_matches_canonical_source(prepare_r_extdata_fixtures):
    for r_name, canonical_rel in R_EXTDATA_COPIES.items():
        r_fixture = R_PACKAGE_DIR / f"inst/extdata/{r_name}"
        canonical = FIXTURES_DIR / canonical_rel
        assert r_fixture.is_file(), "run `make prepare-r-fixtures` to create R extdata copies"
        assert r_fixture.read_bytes() == canonical.read_bytes()


@pytest.mark.r
@skipif_no_rscript
def test_r_extdata_ld_manifest_is_prepared_for_r(prepare_r_extdata_fixtures):
    r_manifest_path = R_PACKAGE_DIR / "inst/extdata/ld/ld_manifest.json"
    r_reference_cache = R_PACKAGE_DIR / "inst/extdata/ld/reference_cache.rds"
    canonical_manifest_path = FIXTURES_DIR / "ld/python/ld_manifest.json"

    assert r_manifest_path.is_file()
    assert r_reference_cache.is_file()

    r_manifest = json.loads(r_manifest_path.read_text())
    canonical_manifest = json.loads(canonical_manifest_path.read_text())
    assert r_manifest.pop("r_reference_cache") == "reference_cache.rds"
    assert len(r_manifest.pop("r_reference_cache_md5")) == 32
    assert r_manifest == canonical_manifest


def test_r_and_python_versions_match():
    import statgen

    assert _r_package_version() == statgen.__version__


def test_r_package_license_matches_repository_license():
    root_license = (REPO_ROOT / "LICENSE").read_text()
    r_license = (R_PACKAGE_DIR / "LICENSE").read_text().splitlines()

    assert "MIT License" in root_license
    assert "Copyright (c) 2026 Oleksandr Frei and contributors" in root_license
    assert "YEAR: 2026" in r_license
    assert "COPYRIGHT HOLDER: Oleksandr Frei and contributors" in r_license


@pytest.mark.r
@skipif_no_rscript
def test_rscript_phase0_source_smoke():
    script = (
        f"source({str(R_PACKAGE_DIR / 'R/verbosity.R')!r}); "
        "stopifnot(identical(get_verbosity(), 'info')); "
        "set_verbosity('quiet'); "
        "stopifnot(identical(get_verbosity(), 'quiet')); "
        "ok <- FALSE; tryCatch(set_verbosity('verbose'), error = function(e) ok <<- TRUE); "
        "stopifnot(ok)"
    )
    result = run_rscript(script)
    assert result.returncode == 0, result.stderr


@pytest.mark.r
@skipif_no_rscript
def test_r_package_install_and_library_smoke(tmp_path, prepare_r_extdata_fixtures):
    deps = run_rscript(
        "pkgs <- c('Matrix', 'jsonlite', 'digest', 'bit64', 'data.table', 'R.utils', 'testthat'); "
        "missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; "
        "cat(paste(missing, collapse = ','))"
    )
    assert deps.returncode == 0, deps.stderr
    missing = deps.stdout.strip()
    if missing:
        pytest.skip(f"R package dependencies unavailable: {missing}")

    lib = tmp_path / "r-lib"
    lib.mkdir()
    install = subprocess.run(
        ["R", "CMD", "INSTALL", "-l", str(lib), str(R_PACKAGE_DIR)],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert install.returncode == 0, install.stderr

    env = os.environ.copy()
    existing_libs = env.get("R_LIBS")
    env["R_LIBS"] = str(lib) if not existing_libs else str(lib) + os.pathsep + existing_libs
    expected_version = _r_package_version()
    result = run_rscript(
        "library(statgen); "
        f"stopifnot(identical(version(), {expected_version!r})); "
        "set_verbosity('quiet'); "
        "stopifnot(identical(get_verbosity(), 'quiet'))",
        env=env,
    )
    assert result.returncode == 0, result.stderr

    testthat_dir = R_PACKAGE_DIR / "tests/testthat"
    result = run_rscript(
        "library(testthat); "
        "library(statgen); "
        f"testthat::test_dir({json.dumps(str(testthat_dir))}, reporter = 'summary')",
        env=env,
        timeout=120,
    )
    assert result.returncode == 0, result.stderr
