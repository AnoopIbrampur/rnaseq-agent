"""Tests for gene lookup statuses and the GSEA ranking/collapse logic."""

from __future__ import annotations

import pandas as pd
import pytest

from pipeline import enrichment, genes

GENE_MAP = pd.DataFrame({
    "gene_id": [1, 2, 3, 4, 5, 6],
    "symbol": ["AAA", "BBB", "CCC", "DDD", "EEE", "AAA"],
    "synonyms": ["A1", "SHARED", "SHARED", None, "ONLYE", None],
    "symbol_source": ["gene_info"] * 5 + ["gene_history:replaced_by_1"],
})
DE = pd.DataFrame({
    "gene_id": [1, 6, 2], "symbol": ["AAA", "AAA", "BBB"],
    "symbol_source": ["gene_info", "gene_history:replaced_by_1", "gene_info"],
    "baseMean": [100.0, 900.0, 50.0], "log2FC": [2.0, 0.1, 0.5], "lfcSE": [0.2, 0.3, 0.1],
    "stat": [10.0, 0.3, 5.0], "padj": [1e-10, 0.9, 1e-4], "significant": [True, False, False],
    "direction": ["up", "ns", "ns"], "note": ["", "", ""],
})
FILTERED = pd.DataFrame({"gene_id": [3], "symbol": ["CCC"], "max_count": [12],
                         "n_samples_passing": [2]})
GSEA = pd.DataFrame({"term": ["SET1", "SET2"], "nes": [2.0, -1.0], "fdr": [0.01, 0.5],
                     "leading_edge": ["AAA;BBB", "BBB"]})
SETS = {"SET1": ["AAA", "BBB"], "SET2": ["BBB"]}


def look(q: str) -> dict:
    return genes.lookup_gene(q, DE, FILTERED, GSEA, SETS, GENE_MAP, fdr_threshold=0.05)


def test_tested_gene_prefers_current_geneid_row() -> None:
    r = look("aaa")
    assert r["status"] == "tested" and r["gene_id"] == 1 and r["direction"] == "up"
    assert r["other_rows_with_same_symbol"] == 1
    assert r["gene_sets"] == [{"term": "SET1", "nes": 2.0, "fdr": 0.01, "significant": True,
                               "in_leading_edge": True}]


def test_padj_significant_but_small_fold_change_is_explained() -> None:
    r = look("BBB")
    assert not r["significant"]
    assert "below the |log2FC|>1 threshold" in r["summary"]


@pytest.mark.parametrize("query,status", [
    ("CCC", "low_count_filtered"), ("DDD", "not_in_count_matrix"),
    ("NOPE", "unknown_symbol"), ("SHARED", "unknown_symbol"),
])
def test_untested_statuses(query: str, status: str) -> None:
    r = look(query)
    assert r["status"] == status
    assert r["summary"]


def test_low_count_message_counts_samples() -> None:
    assert "Only 2 of 8 samples" in look("CCC")["summary"]


def test_alias_resolves_to_official_symbol() -> None:
    r = look("a1")
    assert r["symbol"] == "AAA" and r["resolved_from_alias"] == "a1"
    assert "NCBI" not in r["summary"]


def test_ambiguous_alias_lists_candidates() -> None:
    assert "BBB, CCC" in look("SHARED")["summary"]


def test_collapse_and_rank_are_deterministic() -> None:
    t = pd.DataFrame({
        "symbol": ["X", "X", "Y", "Z", "W"],
        "symbol_source": ["gene_history:discontinued", "gene_info", "gene_info", "gene_info",
                          "gene_info"],
        "baseMean": [999.0, 1.0, 5.0, 5.0, 5.0],
        "stat": [9.0, 1.0, 3.0, 3.0, None],
    })
    r = enrichment.rank_vector(t)
    assert list(r.index) == ["Y", "Z", "X"]  # tie broken by symbol; W (no stat) dropped
    assert r["X"] == 1.0  # current GeneID row wins over higher-baseMean retired row
