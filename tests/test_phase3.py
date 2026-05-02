import json
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

from statgen.annotations import (
    create_annotation,
    create_annotations,
    load_annotations,
    load_annotations_cache,
    save_annotations_cache,
)
from statgen.reference import load_reference
from tests.conftest import FIXTURES_DIR, MATLAB_DIR, run_octave, skipif_no_octave

SHARDED_REF = FIXTURES_DIR / "reference/sharded/@.bim"
ANNO1 = FIXTURES_DIR / "annotations/anno1.bed"
ANNO2 = FIXTURES_DIR / "annotations/anno2.bed"


EXPECTED_MASK = np.array(
    [
        [1, 0],
        [1, 0],
        [1, 0],
        [1, 0],
        [1, 0],
        [0, 1],
        [0, 1],
        [0, 1],
    ],
    dtype=np.uint8,
)


def _write_bed(path: Path, rows: list[tuple[str, int, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for chr_label, start, end in rows:
            f.write(f"{chr_label}\t{start}\t{end}\n")


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def _load_cache_arrays(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        return {k: data[k] for k in data.files}


def test_load_annotations_fixture_masks_and_names():
    reference = load_reference(SHARDED_REF)
    a = load_annotations([ANNO1, ANNO2], reference)

    assert [s.label for s in a.shards] == ["1", "X"]
    assert a.num_snp == 8
    assert a.num_annot == 2
    assert list(a.annonames) == ["anno1", "anno2"]
    assert sparse.isspmatrix_csr(a.annomat)
    np.testing.assert_array_equal(a.annomat.toarray(), EXPECTED_MASK)


def test_load_annotations_accepts_single_path_scalar():
    reference = load_reference(SHARDED_REF)
    a = load_annotations(ANNO1, reference)
    assert list(a.annonames) == ["anno1"]
    np.testing.assert_array_equal(a.annomat.toarray().reshape(-1), EXPECTED_MASK[:, 0])


def test_empty_bed_file_fails(tmp_path):
    reference = load_reference(SHARDED_REF)
    empty = tmp_path / "empty.bed"
    empty.write_text("")
    with pytest.raises(ValueError, match="BED file is empty"):
        load_annotations(empty, reference)


def test_duplicate_bed_basenames_fail(tmp_path):
    reference = load_reference(SHARDED_REF)
    p1 = tmp_path / "set1" / "dup.bed"
    p2 = tmp_path / "set2" / "dup.bed"
    _write_bed(p1, [("1", 99, 100)])
    _write_bed(p2, [("X", 99, 100)])
    with pytest.raises(ValueError, match="duplicate annotation names"):
        load_annotations([p1, p2], reference)


def test_bed_metadata_and_comment_lines_are_ignored(tmp_path):
    reference = load_reference(SHARDED_REF)
    bed = tmp_path / "with_headers.bed"
    _write_text(
        bed,
        "track name=CodingExon description=test\n"
        "browser position chr1:1-1000\n"
        "# comment\n"
        "\n"
        "1\t99\t200\n"
        "1\t299\t400\n",
    )
    a = load_annotations(bed, reference)
    np.testing.assert_array_equal(a.annomat.toarray().reshape(-1), np.array([1, 1, 1, 1, 0, 0, 0, 0], dtype=np.uint8))


def test_metadata_only_bed_file_fails(tmp_path):
    reference = load_reference(SHARDED_REF)
    bed = tmp_path / "metadata_only.bed"
    _write_text(
        bed,
        "track name=foo\n"
        "browser position chr1:1-100\n"
        "# note\n"
        "\n",
    )
    with pytest.raises(ValueError, match="BED file is empty"):
        load_annotations(bed, reference)


def test_whitespace_delimited_bed_rows_supported(tmp_path):
    reference = load_reference(SHARDED_REF)
    bed = tmp_path / "whitespace_delimited.bed"
    _write_text(
        bed,
        "1  99   200\n"
        "1\t299 400\n"
        "X 99\t100\n",
    )
    a = load_annotations(bed, reference)
    np.testing.assert_array_equal(
        a.annomat.toarray().reshape(-1),
        np.array([1, 1, 1, 1, 0, 1, 0, 0], dtype=np.uint8),
    )


def test_boundary_membership_is_start_inclusive_end_exclusive(tmp_path):
    reference = load_reference(SHARDED_REF)
    bed = tmp_path / "boundary.bed"
    _write_bed(
        bed,
        [
            ("1", 99, 100),  # includes bp=100 (pos0=99)
            ("1", 199, 199),  # empty interval, includes nothing
            ("1", 399, 400),  # includes bp=400 (pos0=399)
            ("1", 500, 501),  # excludes bp=500 (pos0=499)
        ],
    )

    a = load_annotations([bed], reference)
    np.testing.assert_array_equal(a.annomat.toarray().reshape(-1), np.array([1, 0, 0, 1, 0, 0, 0, 0], dtype=np.uint8))


def test_adjacent_interval_merge_matches_premerged_intervals(tmp_path):
    reference = load_reference(SHARDED_REF)
    merged = tmp_path / "merged.bed"
    split = tmp_path / "split.bed"
    _write_bed(merged, [("1", 99, 200), ("1", 299, 400)])
    _write_bed(split, [("1", 99, 150), ("1", 150, 200), ("1", 299, 350), ("1", 350, 400)])

    a_merged = load_annotations([merged], reference)
    a_split = load_annotations([split], reference)
    np.testing.assert_array_equal(a_merged.annomat.toarray(), a_split.annomat.toarray())


def test_chromosome_absent_in_bed_yields_zero_mask(tmp_path):
    reference = load_reference(SHARDED_REF)
    bed = tmp_path / "chr2_only.bed"
    _write_bed(bed, [("2", 0, 1000)])

    a = load_annotations([bed], reference)
    np.testing.assert_array_equal(a.annomat.toarray(), np.zeros((reference.num_snp, 1), dtype=np.uint8))


def test_select_annotations_preserves_order_and_rejects_unknown():
    reference = load_reference(SHARDED_REF)
    a = load_annotations([ANNO1, ANNO2], reference)

    sel = a.select_annotations(["anno2", "anno1"])
    assert list(sel.annonames) == ["anno2", "anno1"]
    np.testing.assert_array_equal(sel.annomat.toarray(), EXPECTED_MASK[:, [1, 0]])

    with pytest.raises(ValueError, match="unknown annotation"):
        a.select_annotations(["anno3"])


def test_union_annotations_enforces_compatibility_and_name_collisions():
    reference = load_reference(SHARDED_REF)
    a = load_annotations([ANNO1], reference)
    b = load_annotations([ANNO2], reference)

    u = a.union_annotations(b)
    assert list(u.annonames) == ["anno1", "anno2"]
    np.testing.assert_array_equal(u.annomat.toarray(), EXPECTED_MASK)

    with pytest.raises(ValueError, match="name collision"):
        a.union_annotations(a)

    x_ref = reference.select_shards(["X"])
    x_only = create_annotation(x_ref, np.ones(x_ref.num_snp), "xonly")
    with pytest.raises(ValueError, match="matching shard structure"):
        a.union_annotations(x_only)


def test_create_annotations_and_create_annotation_validation():
    reference = load_reference(SHARDED_REF)

    panel = create_annotations(reference, EXPECTED_MASK, ["anno1", "anno2"])
    np.testing.assert_array_equal(panel.annomat.toarray(), EXPECTED_MASK)

    with pytest.raises(ValueError, match="shape mismatch"):
        create_annotations(reference, EXPECTED_MASK[:, :1], ["anno1", "anno2"])
    with pytest.raises(ValueError, match="non-binary"):
        create_annotations(reference, np.where(EXPECTED_MASK == 1, 2, 0), ["anno1", "anno2"])
    with pytest.raises(ValueError, match="must be unique"):
        create_annotations(reference, EXPECTED_MASK, ["dup", "dup"])

    col = EXPECTED_MASK[:, 0]
    single = create_annotation(reference, col, "anno1")
    np.testing.assert_array_equal(single.annomat.toarray(), col.reshape(-1, 1))

    with pytest.raises(ValueError, match="annovec length mismatch"):
        create_annotation(reference, np.array([1, 0], dtype=np.uint8), "bad")
    with pytest.raises(ValueError, match="must be binary"):
        create_annotation(reference, np.r_[2, np.zeros(reference.num_snp - 1)], "bad")


def test_create_annotation_can_represent_all_snps():
    reference = load_reference(SHARDED_REF)
    panel = create_annotation(reference, np.ones(reference.num_snp, dtype=np.uint8), "all_snps")
    assert list(panel.annonames) == ["all_snps"]
    np.testing.assert_array_equal(panel.annomat.toarray(), np.ones((reference.num_snp, 1), dtype=np.uint8))

    panel_named = create_annotation(reference, np.ones(reference.num_snp, dtype=np.uint8), "custom_all")
    assert list(panel_named.annonames) == ["custom_all"]
    np.testing.assert_array_equal(panel_named.annomat.toarray(), np.ones((reference.num_snp, 1), dtype=np.uint8))


def test_create_annotations_sparse_explicit_zeros_do_not_flip_to_one():
    reference = load_reference(SHARDED_REF)
    # Handcrafted CSR with an explicit zero stored at (0, 0) and a true one at (1, 0).
    # Stored zeros must remain logical zeros after normalization.
    annomat = sparse.csr_matrix(
        (
            np.array([0.0, 1.0], dtype=float),
            np.array([0, 0], dtype=np.int32),
            np.array([0, 1, 2, 2, 2, 2, 2, 2, 2], dtype=np.int32),
        ),
        shape=(reference.num_snp, 1),
    )

    panel = create_annotations(reference, annomat, ["anno"])
    dense = panel.annomat.toarray().reshape(-1)
    assert dense[0] == 0
    assert dense[1] == 1
    assert np.count_nonzero(dense) == 1


def test_cache_roundtrip_subset_and_compatibility(tmp_path):
    reference = load_reference(SHARDED_REF)
    a = load_annotations([ANNO1, ANNO2], reference)
    cache = tmp_path / "annotations_cache.npz"

    save_annotations_cache(a, cache)
    loaded = load_annotations_cache(cache)

    assert sparse.isspmatrix_csr(loaded.annomat)
    np.testing.assert_array_equal(loaded.annomat.toarray(), EXPECTED_MASK)
    assert reference.is_object_compatible(loaded) is True

    x_only = load_annotations_cache(cache, shards=["X"])
    assert [sh.label for sh in x_only.shards] == ["X"]
    np.testing.assert_array_equal(x_only.annomat.toarray(), EXPECTED_MASK[5:, :])


def test_cache_validation_and_post_load_compatibility_check(tmp_path):
    reference = load_reference(SHARDED_REF)
    a = load_annotations([ANNO1, ANNO2], reference)
    cache = tmp_path / "annotations_cache.npz"
    save_annotations_cache(a, cache)

    arrays = _load_cache_arrays(cache)
    meta = json.loads(bytes(arrays["_meta"]).decode())

    meta["schema"] = "annotations_cache/bad"
    arrays["_meta"] = np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8)
    bad_schema = tmp_path / "bad_schema.npz"
    np.savez_compressed(bad_schema, **arrays)
    with pytest.raises(ValueError, match="Unsupported annotations cache schema"):
        load_annotations_cache(bad_schema)

    arrays = _load_cache_arrays(cache)
    meta = json.loads(bytes(arrays["_meta"]).decode())
    meta["shard_checksums"] = meta["shard_checksums"][:-1]
    arrays["_meta"] = np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8)
    bad_len = tmp_path / "bad_len.npz"
    np.savez_compressed(bad_len, **arrays)
    with pytest.raises(ValueError, match="shard_labels and shard_checksums length mismatch"):
        load_annotations_cache(bad_len)

    arrays = _load_cache_arrays(cache)
    meta = json.loads(bytes(arrays["_meta"]).decode())
    meta["shard_checksums"][0] = "deadbeef" * 4
    arrays["_meta"] = np.frombuffer(json.dumps(meta).encode(), dtype=np.uint8)
    bad_chk = tmp_path / "bad_chk.npz"
    np.savez_compressed(bad_chk, **arrays)
    loaded_bad = load_annotations_cache(bad_chk)
    assert reference.is_object_compatible(loaded_bad) is False


def test_cross_language_cache_not_supported(tmp_path):
    reference = load_reference(SHARDED_REF)
    a = load_annotations([ANNO1, ANNO2], reference)

    py_cache = tmp_path / "annotations_cache_py.npz"
    save_annotations_cache(a, py_cache)
    octave_result = run_octave(
        _octave_script(
            f"try; statgen.load_annotations_cache('{py_cache}'); fprintf('NOFAIL\\n'); catch; fprintf('FAIL\\n'); end"
        )
    )
    assert octave_result.returncode == 0, octave_result.stderr
    assert octave_result.stdout.strip() == "FAIL"

    mat_cache = tmp_path / "annotations_cache_mat.mat"
    octave_result = run_octave(
        _octave_script(
            f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
            f"a = statgen.load_annotations({{[fixture_dir '/annotations/anno1.bed'], [fixture_dir '/annotations/anno2.bed']}}, ref); "
            f"statgen.save_annotations_cache(a, '{mat_cache}'); "
            "fprintf('OK\\n');"
        )
    )
    assert octave_result.returncode == 0, octave_result.stderr
    assert octave_result.stdout.strip() == "OK"
    with pytest.raises(Exception):
        load_annotations_cache(mat_cache)


def _octave_script(expr: str) -> str:
    matlab_dir = str(MATLAB_DIR)
    fixture_dir = str(FIXTURES_DIR)
    return f"addpath('{matlab_dir}'); fixture_dir = '{fixture_dir}'; " + expr


@pytest.mark.octave
@skipif_no_octave
def test_octave_load_annotations_fixture_masks_and_names():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "a = statgen.load_annotations({[fixture_dir '/annotations/anno1.bed'], [fixture_dir '/annotations/anno2.bed']}, ref); "
        "fprintf('%d\\n', a.num_snp); "
        "fprintf('%d\\n', a.num_annot); "
        "fprintf('%s,%s\\n', a.annonames{1}, a.annonames{2}); "
        "fprintf('%d,%d\\n', nnz(a.annomat(:,1)), nnz(a.annomat(:,2))); "
        "M = full(a.annomat); "
        "fprintf('%d', M(1,1)); fprintf('%d', M(5,1)); fprintf('%d', M(6,1)); fprintf('%d', M(6,2)); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "8"
    assert lines[1] == "2"
    assert lines[2] == "anno1,anno2"
    assert lines[3] == "5,3"
    assert lines[4] == "1101"


@pytest.mark.octave
@skipif_no_octave
def test_octave_annotations_cache_roundtrip_and_subset(tmp_path):
    cache = tmp_path / "annotations_cache.mat"
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "a = statgen.load_annotations({[fixture_dir '/annotations/anno1.bed'], [fixture_dir '/annotations/anno2.bed']}, ref); "
        f"statgen.save_annotations_cache(a, '{cache}'); "
        f"b = statgen.load_annotations_cache('{cache}'); "
        f"x = statgen.load_annotations_cache('{cache}', {{'X'}}); "
        "fprintf('%d\\n', ref.is_object_compatible(b)); "
        "fprintf('%d\\n', b.num_snp); "
        "fprintf('%d\\n', x.num_snp); "
        "fprintf('%s\\n', x.shards{1}.label);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "1"
    assert lines[1] == "8"
    assert lines[2] == "3"
    assert lines[3] == "X"


@pytest.mark.octave
@skipif_no_octave
def test_octave_empty_bed_fails(tmp_path):
    empty = tmp_path / "empty.bed"
    empty.write_text("")
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"try; statgen.load_annotations('{empty}', ref); fprintf('NOFAIL\\n'); catch; fprintf('FAIL\\n'); end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "FAIL"


@pytest.mark.octave
@skipif_no_octave
def test_octave_create_annotation_can_represent_all_snps():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "a = statgen.create_annotation(ref, ones(ref.num_snp,1), 'all_snps'); "
        "b = statgen.create_annotation(ref, ones(ref.num_snp,1), 'custom_all'); "
        "fprintf('%s\\n', a.annonames{1}); "
        "fprintf('%d\\n', nnz(a.annomat)); "
        "fprintf('%s\\n', b.annonames{1}); "
        "fprintf('%d\\n', nnz(b.annomat));"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "all_snps"
    assert lines[1] == "8"
    assert lines[2] == "custom_all"
    assert lines[3] == "8"


@pytest.mark.octave
@skipif_no_octave
def test_octave_duplicate_bed_basenames_fail(tmp_path):
    p1 = tmp_path / "set1" / "dup.bed"
    p2 = tmp_path / "set2" / "dup.bed"
    _write_bed(p1, [("1", 99, 100)])
    _write_bed(p2, [("X", 99, 100)])

    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"try; statgen.load_annotations({{'{p1}','{p2}'}}, ref); fprintf('NOFAIL\\n'); catch; fprintf('FAIL\\n'); end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "FAIL"


@pytest.mark.octave
@skipif_no_octave
def test_octave_bed_metadata_and_comment_lines_are_ignored(tmp_path):
    bed = tmp_path / "with_headers_octave.bed"
    _write_text(
        bed,
        "track name=CodingExon description=test\n"
        "browser position chr1:1-1000\n"
        "# comment\n"
        "\n"
        "1\t99\t200\n"
        "1\t299\t400\n",
    )
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"a = statgen.load_annotations('{bed}', ref); "
        "v = full(a.annomat(:,1)); "
        "fprintf('%d', v(1)); fprintf('%d', v(2)); fprintf('%d', v(3)); fprintf('%d', v(4)); "
        "fprintf('%d', v(5)); fprintf('%d', v(6)); fprintf('%d', v(7)); fprintf('%d', v(8)); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "11110000"


@pytest.mark.octave
@skipif_no_octave
def test_octave_metadata_only_bed_file_fails(tmp_path):
    bed = tmp_path / "metadata_only_octave.bed"
    _write_text(
        bed,
        "track name=foo\n"
        "browser position chr1:1-100\n"
        "# note\n"
        "\n",
    )
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"try; statgen.load_annotations('{bed}', ref); fprintf('NOFAIL\\n'); catch; fprintf('FAIL\\n'); end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "FAIL"


@pytest.mark.octave
@skipif_no_octave
def test_octave_whitespace_delimited_bed_rows_supported(tmp_path):
    bed = tmp_path / "whitespace_delimited_octave.bed"
    _write_text(
        bed,
        "1  99   200\n"
        "1\t299 400\n"
        "X 99\t100\n",
    )
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"a = statgen.load_annotations('{bed}', ref); "
        "v = full(a.annomat(:,1)); "
        "fprintf('%d', v(1)); fprintf('%d', v(2)); fprintf('%d', v(3)); fprintf('%d', v(4)); "
        "fprintf('%d', v(5)); fprintf('%d', v(6)); fprintf('%d', v(7)); fprintf('%d', v(8)); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "11110100"


@pytest.mark.octave
@skipif_no_octave
def test_octave_boundary_membership_start_inclusive_end_exclusive(tmp_path):
    bed = tmp_path / "boundary_octave.bed"
    _write_bed(
        bed,
        [
            ("1", 99, 100),
            ("1", 199, 199),
            ("1", 399, 400),
            ("1", 500, 501),
        ],
    )
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"a = statgen.load_annotations('{bed}', ref); "
        "v = full(a.annomat(:,1)); "
        "fprintf('%d', v(1)); fprintf('%d', v(2)); fprintf('%d', v(3)); fprintf('%d', v(4)); "
        "fprintf('%d', v(5)); fprintf('%d', v(6)); fprintf('%d', v(7)); fprintf('%d', v(8)); fprintf('\\n');"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "10010000"


@pytest.mark.octave
@skipif_no_octave
def test_octave_adjacent_interval_merge_matches_premerged(tmp_path):
    merged = tmp_path / "merged_octave.bed"
    split = tmp_path / "split_octave.bed"
    _write_bed(merged, [("1", 99, 200), ("1", 299, 400)])
    _write_bed(split, [("1", 99, 150), ("1", 150, 200), ("1", 299, 350), ("1", 350, 400)])
    script = _octave_script(
        f"ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        f"a = statgen.load_annotations('{merged}', ref); "
        f"b = statgen.load_annotations('{split}', ref); "
        "eq = isequal(full(a.annomat), full(b.annomat)); "
        "fprintf('%d\\n', eq);"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "1"


@pytest.mark.octave
@skipif_no_octave
def test_octave_union_annotations_happy_path_and_name_collision():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "a = statgen.load_annotations([fixture_dir '/annotations/anno1.bed'], ref); "
        "b = statgen.load_annotations([fixture_dir '/annotations/anno2.bed'], ref); "
        "u = a.union_annotations(b); "
        "fprintf('%s,%s\\n', u.annonames{1}, u.annonames{2}); "
        "fprintf('%d,%d\\n', nnz(u.annomat(:,1)), nnz(u.annomat(:,2))); "
        "try; a.union_annotations(a); fprintf('NOFAIL\\n'); catch; fprintf('FAIL\\n'); end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "anno1,anno2"
    assert lines[1] == "5,3"
    assert lines[2] == "FAIL"


@pytest.mark.octave
@skipif_no_octave
def test_octave_create_annotations_and_create_annotation_validation():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "A = [ones(ref.num_snp,1), zeros(ref.num_snp,1)]; "
        "ok = statgen.create_annotations(ref, A, {'a','b'}); "
        "fprintf('%d\\n', ok.num_annot); "
        "try; statgen.create_annotations(ref, ones(ref.num_snp,1), {'a','b'}); fprintf('NOFAIL1\\n'); catch; fprintf('FAIL1\\n'); end; "
        "try; statgen.create_annotations(ref, 2*ones(ref.num_snp,2), {'a','b'}); fprintf('NOFAIL2\\n'); catch; fprintf('FAIL2\\n'); end; "
        "try; statgen.create_annotations(ref, A, {'dup','dup'}); fprintf('NOFAIL3\\n'); catch; fprintf('FAIL3\\n'); end; "
        "single = statgen.create_annotation(ref, ones(ref.num_snp,1), 'single'); fprintf('%d\\n', single.num_annot); "
        "try; statgen.create_annotation(ref, ones(2,1), 'bad'); fprintf('NOFAIL4\\n'); catch; fprintf('FAIL4\\n'); end; "
        "v = ones(ref.num_snp,1); v(1)=2; try; statgen.create_annotation(ref, v, 'bad'); fprintf('NOFAIL5\\n'); catch; fprintf('FAIL5\\n'); end;"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "2"
    assert lines[1] == "FAIL1"
    assert lines[2] == "FAIL2"
    assert lines[3] == "FAIL3"
    assert lines[4] == "1"
    assert lines[5] == "FAIL4"
    assert lines[6] == "FAIL5"


@pytest.mark.octave
@skipif_no_octave
def test_octave_select_annotations_and_unknown_name_error():
    script = _octave_script(
        "ref = statgen.load_reference([fixture_dir '/reference/sharded/@.bim']); "
        "a = statgen.load_annotations({[fixture_dir '/annotations/anno1.bed'], [fixture_dir '/annotations/anno2.bed']}, ref); "
        "b = a.select_annotations({'anno2','anno1'}); "
        "fprintf('%s,%s\\n', b.annonames{1}, b.annonames{2}); "
        "try; a.select_annotations({'missing'}); fprintf('NOFAIL\\n'); catch; fprintf('FAIL\\n'); end"
    )
    result = run_octave(script)
    assert result.returncode == 0, result.stderr
    lines = result.stdout.strip().splitlines()
    assert lines[0] == "anno2,anno1"
    assert lines[1] == "FAIL"
