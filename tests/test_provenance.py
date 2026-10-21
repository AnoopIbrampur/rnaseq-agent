"""Tests for run records, cache keys, cache hits, and reproduce."""

from __future__ import annotations

import json
import shutil

import pytest

from pipeline import cache, config, provenance

QC_PARAMS = {"qc_rules": config.QC_RULES, "low_count": config.LOW_COUNT_FILTER}
REQUIRED_FIELDS = {
    "run_id", "tool", "timestamp_utc", "inputs", "parameters", "software",
    "code_sha256", "git", "outputs", "output_sha256", "cache_hit", "cache_key",
    "duration_s", "upstream_runs", "summary", "result",
}


@pytest.fixture
def dirs(tmp_path):
    return {"runs_dir": tmp_path / "runs", "results_root": tmp_path / "results"}


@pytest.fixture
def qc_inputs():
    return {"counts": config.COUNTS_FILE, "metadata": config.SERIES_MATRIX_FILE}


def test_cache_key_is_stable_and_sensitive() -> None:
    base = dict(tool="t", input_checksums={"a": "1"}, parameters={"x": 1, "y": [1, 2]},
                software={"p": "1.0"}, code_sha="c")
    k = cache.cache_key(**base)
    assert k == cache.cache_key(**{**base, "parameters": {"y": [1, 2], "x": 1}})  # key order
    for change in ({"input_checksums": {"a": "2"}}, {"parameters": {"x": 2, "y": [1, 2]}},
                   {"software": {"p": "1.1"}}, {"code_sha": "d"}, {"tool": "u"}):
        assert cache.cache_key(**{**base, **change}) != k, change


def test_record_is_complete_and_on_disk(dirs, qc_inputs) -> None:
    rec = provenance.execute("run_qc", qc_inputs, QC_PARAMS, **dirs)
    assert REQUIRED_FIELDS <= set(rec)
    on_disk = json.loads((dirs["runs_dir"] / f"{rec['run_id']}.json").read_text())
    assert on_disk == rec
    for name, inp in rec["inputs"].items():
        assert len(inp["sha256"]) == 64 and inp["path"]
    for name, path in rec["outputs"].items():
        assert provenance.absolute(path).exists()
        assert len(rec["output_sha256"][name]) == 64
    assert {"python", "pydeseq2", "gseapy", "pandas"} <= set(rec["software"])
    assert rec["result"]["passed"] is True


def test_second_call_is_cache_hit_with_new_run_id(dirs, qc_inputs) -> None:
    a = provenance.execute("run_qc", qc_inputs, QC_PARAMS, **dirs)
    b = provenance.execute("run_qc", qc_inputs, QC_PARAMS, **dirs)
    assert not a["cache_hit"] and b["cache_hit"]
    assert a["run_id"] != b["run_id"]
    assert a["cache_key"] == b["cache_key"] and a["result"] == b["result"]


def test_changed_parameter_misses_cache(dirs, qc_inputs) -> None:
    a = provenance.execute("run_qc", qc_inputs, QC_PARAMS, **dirs)
    params = {**QC_PARAMS, "low_count": {"min_count": 5, "min_samples": 4}}
    b = provenance.execute("run_qc", qc_inputs, params, **dirs)
    assert not b["cache_hit"] and a["cache_key"] != b["cache_key"]


def test_reproduce_qc_succeeds(dirs, qc_inputs) -> None:
    rec = provenance.execute("run_qc", qc_inputs, QC_PARAMS, **dirs)
    rep = provenance.reproduce(rec["run_id"], runs_dir=dirs["runs_dir"])
    assert rep["ok"], rep
    assert (dirs["runs_dir"] / f"{rec['run_id']}.reproduce.json").exists()


def test_reproduce_detects_changed_input(dirs, tmp_path) -> None:
    meta = tmp_path / "in" / config.SERIES_MATRIX_FILE.name
    meta.parent.mkdir()
    shutil.copy(config.SERIES_MATRIX_FILE, meta)
    rec = provenance.execute("run_qc", {"counts": config.COUNTS_FILE, "metadata": meta},
                             QC_PARAMS, **dirs)
    with open(meta, "ab") as f:
        f.write(b"\n")
    rep = provenance.reproduce(rec["run_id"], runs_dir=dirs["runs_dir"])
    assert not rep["ok"] and "changed since the run" in rep["error"]


def test_reproduce_detects_tampered_output(dirs, qc_inputs) -> None:
    rec = provenance.execute("run_qc", qc_inputs, QC_PARAMS, **dirs)
    path = dirs["runs_dir"] / f"{rec['run_id']}.json"
    tampered = json.loads(path.read_text())
    tampered["result"]["metrics"]["genes_passing_filter"] += 1
    path.write_text(json.dumps(tampered))
    rep = provenance.reproduce(rec["run_id"], runs_dir=dirs["runs_dir"])
    assert not rep["ok"]
    assert any("genes_passing_filter" in d for c in rep["checks"] for d in c["detail"])


def test_unknown_run_id_has_helpful_error(dirs) -> None:
    dirs["runs_dir"].mkdir(parents=True)
    with pytest.raises(FileNotFoundError, match="Known run ids"):
        provenance.load_record("nope", dirs["runs_dir"])


def test_values_match_tolerance() -> None:
    assert provenance.values_match({"a": 1.0, "b": [1, "x"]}, {"a": 1.0 + 1e-9, "b": [1, "x"]}) == []
    assert provenance.values_match({"a": 1.0}, {"a": 1.01})
    assert provenance.values_match({"a": 1}, {"a": 1, "b": 2})
