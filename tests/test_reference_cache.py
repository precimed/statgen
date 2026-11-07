"""Python reference-cache schema migration and string-storage regressions."""
import json

import numpy as np
import pytest

from statgen import reference as ref_mod
from statgen.reference import ReferencePanel, ReferenceShard, load_reference, load_reference_cache
from tests.conftest import FIXTURES_DIR


FIELDS = ('chr', 'snp', 'bp', 'a1', 'a2', 'a1_hash64', 'a2_hash64')


def assert_same(actual, expected):
    assert [s.label for s in actual.shards] == [s.label for s in expected.shards]
    for left, right in zip(actual.shards, expected.shards):
        assert left.checksum == right.checksum
        for field in FIELDS:
            np.testing.assert_array_equal(getattr(left, field), getattr(right, field))
        for field in ('chr', 'snp', 'a1', 'a2'):
            assert getattr(left, field).dtype == object


def cache_arrays(tmp_path):
    panel = load_reference(FIXTURES_DIR / 'reference/sharded/@.bim')
    path = tmp_path / 'ref.npz'
    panel.save_cache(path)
    with np.load(path, allow_pickle=False) as data:
        arrays = {key: data[key] for key in data.files}
    return path, arrays


@pytest.mark.parametrize('schema', ['0_1', '0_2'])
def test_schema_fixtures_match_sources(schema):
    cached = load_reference_cache(FIXTURES_DIR / f'cache/reference_python_reference_cache_{schema}.npz')
    source = load_reference(FIXTURES_DIR / 'reference/sharded/@.bim')
    assert_same(cached, source)
    selected = load_reference_cache(
        FIXTURES_DIR / f'cache/reference_python_reference_cache_{schema}.npz', shards=['X'])
    assert_same(selected, source.select_shards(['X']))


@pytest.mark.parametrize('fixture', [
    'cache/reference_python_reference_cache_0_1.npz',
    'long_alleles/legacy_python.npz',
])
def test_migration_preserves_stored_hashes_without_recomputing(tmp_path, monkeypatch, fixture):
    old = load_reference_cache(FIXTURES_DIR / fixture)
    def fail(*args, **kwargs):
        raise AssertionError('cache migration must preserve stored identities')
    monkeypatch.setattr(ref_mod, 'allele_hash64', fail)
    monkeypatch.setattr(ref_mod, '_checksum_from_arrays', fail)
    path = tmp_path / 'new.npz'
    old.save_cache(path)
    with np.load(path, allow_pickle=False) as data:
        assert json.loads(data['_meta'].tobytes())['schema'] == 'reference_cache/0.2'
        assert 's0_a1' not in data
    assert_same(load_reference_cache(path), old)


def test_empty_shard_and_selection(tmp_path):
    empty = ReferenceShard('1', [], [], [], [], [])
    populated = ReferenceShard('X', ['X'], ['rsé'], [10], ['A'], ['G'])
    panel = ReferencePanel([empty, populated])
    path = tmp_path / 'empty.npz'
    panel.save_cache(path)
    assert_same(load_reference_cache(path), panel)
    with np.load(path) as data:
        assert data['s0_snp_text'].size == data['s0_a1_text'].size == data['s0_a2_text'].size == 0
        arrays = {key: data[key] for key in data.files}
    # Selection must not decode unselected payloads.
    arrays['s0_a1_text'] = np.array([255], dtype=np.uint8)
    np.savez(path, **arrays)
    selected = load_reference_cache(path, shards=['X'])
    assert_same(selected, panel.select_shards(['X']))
    assert selected.select_shards(['X']).shards[0] is selected.shards[0]


@pytest.mark.parametrize('field,value,error', [
    ('s0_a1_text', np.array([255], dtype=np.uint8), 'UTF-8'),
    ('s0_a1_text', np.array([], dtype=np.uint8), 'length mismatch'),
    ('s0_a1_text', np.frombuffer(b'A\nC\nG\nT\nA\n', dtype=np.uint8), 'length mismatch'),
    ('s0_snp_text', np.array([65], dtype=np.int64), 'uint8'),
    ('s0_a2_text', np.array([[65]], dtype=np.uint8), 'one-dimensional'),
    ('s0_bp', np.arange(5, dtype=np.int32), 'int64'),
    ('s0_bp', np.zeros((5, 1), dtype=np.int64), 'one-dimensional'),
    ('s0_a1_hash64', np.zeros(5, dtype=np.int64), 'uint64'),
    ('s0_a2_hash64', np.zeros(4, dtype=np.uint64), 'lengths mismatch'),
    ('s0_chr', np.array(['1']), 'lengths mismatch'),
    ('s0_snp_text', None, 'missing s0_snp_text'),
    ('_meta', np.array([1], dtype=np.int64), 'uint8'),
])
def test_malformed_payload(tmp_path, field, value, error):
    path, arrays = cache_arrays(tmp_path)
    if value is None:
        del arrays[field]
    else:
        arrays[field] = value
    np.savez(path, **arrays)
    with pytest.raises(ValueError, match=error):
        load_reference_cache(path)


@pytest.mark.parametrize('metadata', [
    {'schema': 'reference_cache/0.2'},
    {'schema': 'reference_cache/0.2', 'shard_labels': ['1', '1'], 'shard_checksums': ['x', 'x']},
    {'schema': 'reference_cache/0.2', 'shard_labels': ['1'], 'shard_checksums': []},
])
def test_required_metadata(tmp_path, metadata):
    path = tmp_path / 'bad.npz'
    np.savez(path, _meta=np.frombuffer(json.dumps(metadata).encode(), dtype=np.uint8))
    with pytest.raises(ValueError, match='Invalid reference cache'):
        load_reference_cache(path)


@pytest.mark.parametrize('schema', ['reference_cache/0.3', 'reference_cache/99', 'other/0.2'])
def test_unsupported_schema_before_payload(tmp_path, schema):
    path = tmp_path / 'future.npz'
    np.savez(path, _meta=np.frombuffer(json.dumps({'schema': schema}).encode(), dtype=np.uint8))
    with pytest.raises(ValueError, match=f'Unsupported reference cache schema: {schema!r}'):
        load_reference_cache(path)


def test_long_outlier_payload_sizes_roundtrip_and_determinism(tmp_path):
    count, length = 2000, 193000
    snp = ['rs'] * count
    a1 = ['A'] * count
    snp[-1] = 'id_é' + 'x' * length
    a1[-1] = 'C' * length
    shard = ReferenceShard('1', ['1'] * count, snp, np.arange(count), a1, ['G'] * count)
    other = ReferenceShard('X', ['X'], ['other'], [10], ['A'], ['G'])
    panel = ReferencePanel([shard, other])
    path = tmp_path / 'long.npz'
    panel.save_cache(path)
    assert_same(load_reference_cache(path), panel)
    with np.load(path, allow_pickle=False) as data:
        assert data['s0_a1_text'].nbytes == length + 2 * (count - 1)
        assert data['s0_snp_text'].nbytes == len(('\n'.join(snp)).encode())
    again = tmp_path / 'again.npz'
    panel.save_cache(again)
    assert path.read_bytes() == again.read_bytes()


@pytest.mark.parametrize('byteorder', ['<', '>'])
def test_numeric_vectors_accept_either_byte_order(tmp_path, byteorder):
    path, arrays = cache_arrays(tmp_path)
    expected = load_reference_cache(path)
    for field in ('bp', 'a1_hash64', 'a2_hash64'):
        key = 's0_' + field
        arrays[key] = arrays[key].astype(arrays[key].dtype.newbyteorder(byteorder))
    np.savez(path, **arrays)
    assert_same(load_reference_cache(path), expected)
