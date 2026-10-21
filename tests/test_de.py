"""Tests for the DE tool: direction of effect, expected genes, output schema."""

from __future__ import annotations

import pandas as pd
import pytest

from pipeline import data, de

EXPECTED_UP = ["FKBP5", "TSC22D3", "PER1", "KLF15", "DUSP1", "ZBTB16", "CRISPLD2"]


@pytest.fixture(scope="module")
def de_run(tmp_path_factory):
    ds = data.load_dataset()
    out = tmp_path_factory.mktemp("de")
    return de.run_de(ds.counts, ds.metadata, ds.gene_map, out_dir=out)


def test_reference_level_is_untreated() -> None:
    meta = pd.DataFrame({"cell": ["a", "a"], "dex": ["trt", "untrt"]})
    fixed = de._design_metadata(meta)
    assert list(fixed["dex"].cat.categories) == ["untrt", "trt"]


def test_known_dex_genes_are_up(de_run) -> None:
    t = pd.read_csv(de_run["artifacts"]["table"]).set_index("symbol")
    up = [g for g in EXPECTED_UP if t.loc[g, "direction"] == "up"]
    # Phase 2 acceptance criterion: at least 5 of the literature genes.
    assert len(up) >= 5, f"only {up} up"


def test_output_schema_and_counts_consistent(de_run) -> None:
    t = pd.read_csv(de_run["artifacts"]["table"])
    for col in ["gene_id", "symbol", "baseMean", "log2FC", "padj", "stat", "direction", "note"]:
        assert col in t.columns
    assert de_run["n_sig_up"] == (t["direction"] == "up").sum()
    assert de_run["n_sig_down"] == (t["direction"] == "down").sum()
    sig = t[t["significant"]]
    assert (sig["padj"] < 0.05).all() and (sig["log2FC"].abs() > 1).all()
    assert de_run["n_low_count_filtered"] + de_run["n_genes_modelled"] == 39376
