"""Tests for QC gate logic, on synthetic inputs and on real/corrupted data."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pipeline import config, data, qc


@pytest.fixture(scope="module")
def ds() -> data.Dataset:
    return data.load_dataset()


def test_library_size_rule_flags_only_small_samples() -> None:
    lib = pd.Series({"a": 100, "b": 100, "c": 100, "d": 49})
    flags = qc.library_size_flags(lib, 0.5)
    assert [f["sample"] for f in flags] == ["d"]
    assert flags[0]["rule"] == "min_library_fraction"
    assert flags[0]["threshold"] == 50.0


def test_library_size_rule_threshold_is_strict() -> None:
    lib = pd.Series({"a": 100, "b": 100, "c": 50})
    assert qc.library_size_flags(lib, 0.5) == []


def test_pca_outlier_rule_flags_far_sample() -> None:
    scores = pd.DataFrame(
        {"PC1": [0.0, 1.0, 0.5, 0.2, 10.0, 11.0, 10.5, 40.0],
         "PC2": [0.0, 0.5, 1.0, 0.3, 0.0, 0.5, 1.0, 30.0]},
        index=list("abcdefgh"),
    )
    groups = pd.Series(["x"] * 4 + ["y"] * 4, index=scores.index)
    flags, _ = qc.pca_outlier_flags(scores, groups, 3.0)
    assert [f["sample"] for f in flags] == ["h"]


def test_pca_outlier_rule_no_flags_when_spread_is_zero() -> None:
    scores = pd.DataFrame({"PC1": [1.0, 1.0, 2.0, 2.0], "PC2": [0.0] * 4}, index=list("abcd"))
    groups = pd.Series(["x", "x", "y", "y"], index=scores.index)
    flags, _ = qc.pca_outlier_flags(scores, groups, 3.0)
    assert flags == []


def test_eta_squared_extremes() -> None:
    g = pd.Series(["a", "a", "b", "b"])
    assert qc.eta_squared(pd.Series([0.0, 0.0, 1.0, 1.0]), g) == pytest.approx(1.0)
    assert qc.eta_squared(pd.Series([0.0, 1.0, 0.0, 1.0]), g) == pytest.approx(0.0)


def test_low_count_mask() -> None:
    counts = pd.DataFrame({"s1": [10, 9, 100], "s2": [10, 100, 0], "s3": [0, 0, 0]})
    assert qc.low_count_mask(counts, 10, 2).tolist() == [True, False, False]


def test_real_data_passes_and_separates_by_both(ds, tmp_path) -> None:
    r = qc.run_qc(ds.counts, ds.metadata, out_dir=tmp_path)
    assert r["passed"], r["summary"]
    assert r["flags"] == []
    sep = r["metrics"]["separation"]
    assert sep["dex_pcs"] == ["PC1"] and sep["cell_pcs"] == ["PC2"]
    for path in r["artifacts"].values():
        assert (tmp_path / path.split("/")[-1]).exists()


@pytest.mark.parametrize("kind,rule", [
    ("low_library", "min_library_fraction"),
    ("scrambled_sample", "pca_outlier_mad"),
])
def test_corrupted_data_fails_gate(ds, tmp_path, kind, rule) -> None:
    sample = ds.counts.columns[3]
    bad = data.corrupt_counts(ds.counts, kind, sample)
    r = qc.run_qc(bad, ds.metadata, out_dir=tmp_path)
    assert not r["passed"]
    assert [(f["sample"], f["rule"]) for f in r["flags"]] == [(sample, rule)]
    assert "FAILED" in r["summary"]


def test_corrupt_counts_is_deterministic_and_rejects_unknown(ds) -> None:
    s = ds.counts.columns[0]
    a = data.corrupt_counts(ds.counts, "low_library", s)
    b = data.corrupt_counts(ds.counts, "low_library", s)
    assert a.equals(b)
    assert np.array_equal(a.drop(columns=s), ds.counts.drop(columns=s))
    with pytest.raises(ValueError, match="Unknown corruption"):
        data.corrupt_counts(ds.counts, "bogus", s)


def test_rules_have_explanations() -> None:
    for key in ("min_library_fraction", "pca_outlier_mad", "separation_r2"):
        assert config.QC_RULES[f"{key}_why"]
