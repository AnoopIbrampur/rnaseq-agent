"""Tests for data loading, checksum verification, and gene mapping."""

from __future__ import annotations

import shutil

import pytest

from pipeline import config, data


@pytest.fixture(scope="module")
def ds() -> data.Dataset:
    return data.load_dataset()


def test_raw_checksums_match_config() -> None:
    for name in config.EXPECTED_SHA256:
        data.verify_raw_file(config.DATA_RAW / name)


def test_checksum_mismatch_raises(tmp_path) -> None:
    bad = tmp_path / config.SERIES_MATRIX_FILE.name
    shutil.copy(config.SERIES_MATRIX_FILE, bad)
    with open(bad, "ab") as f:
        f.write(b"tampered")
    with pytest.raises(data.DataIntegrityError, match="Checksum mismatch"):
        data.verify_raw_file(bad)


def test_missing_file_raises(tmp_path) -> None:
    with pytest.raises(data.DataIntegrityError, match="Missing raw file"):
        data.verify_raw_file(tmp_path / "nope.tsv.gz")


def test_metadata_is_paired_design(ds: data.Dataset) -> None:
    m = ds.metadata
    assert len(m) == 8
    assert list(m["dex"].cat.categories) == ["untrt", "trt"]
    # Every donor cell line contributes exactly one untreated and one treated sample.
    assert (m.groupby("cell", observed=True)["dex"].value_counts().unstack() == 1).all().all()


def test_counts_align_with_metadata(ds: data.Dataset) -> None:
    assert list(ds.counts.columns) == list(ds.metadata.index)
    assert ds.counts.shape == (39376, 8)
    assert (ds.counts >= 0).all().all()


def test_every_count_geneid_has_a_symbol(ds: data.Dataset) -> None:
    gm = ds.gene_map.set_index("gene_id")
    assert not ds.gene_map["gene_id"].duplicated().any()
    assert ds.counts.index.isin(gm.index).all()
    assert gm.loc[2289, "symbol"] == "FKBP5"
