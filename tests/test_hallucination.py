"""Tests for the claim extractor and grounding checks."""

from __future__ import annotations

import pytest

from eval import hallucination as H
from eval.run_eval import claimed_status


def test_split_context() -> None:
    f, c = H.split_context("FKBP5 is up [run:abcd1234].\n\nContext (not computed here): GC target.")
    assert "FKBP5" in f and "GC target" in c and "Context" not in f


def test_extract_genes_skips_stopwords_sets_and_citations() -> None:
    text = ("FKBP5 and DUSP1 are higher; PCA and QC passed; NES +1.98 for "
            "HALLMARK_IL2_STAT5_SIGNALING [run:abcd1234]. GC response.")
    assert H.extract_genes(text) == ["FKBP5", "DUSP1"]


def test_extract_pathways_with_and_without_prefix() -> None:
    assert H.extract_pathways("HALLMARK_ADIPOGENESIS and TNFA_SIGNALING_VIA_NFKB") == [
        "HALLMARK_ADIPOGENESIS", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"]
    assert H.extract_pathways("NOT_A_REAL_SET") == []


def test_extract_numbers_ignores_symbols_lists_and_log2() -> None:
    text = ("1. KLF15: log2FC +4.39, padj 1.18e-113 [run:abcd1234]\n"
            "2. 458 genes up in GSE52778, PC1 43% of variance, 3.2 × 10^-5")
    assert H.extract_numbers(text) == ["+4.39", "1.18e-113", "458", "43", "3.2e-5"]


@pytest.mark.parametrize("written,sources,ok", [
    ("3.96", [3.962217], True),
    ("4.0", [3.962217], True),
    ("4", [3.962217], True),          # integer written, half-unit tolerance
    ("3.9", [3.962217], False),
    ("2.13e-26", [2.132311e-26], True),
    ("2.1e-26", [2.132311e-26], True),
    ("5e-26", [2.132311e-26], False),
    ("-3.65", [-3.646868], True),
    ("3.65", [-3.646868], True),      # magnitude: "lower by 3.65"
    ("458", [458.0], True),
    ("459", [458.0], False),
    ("1,234", [1234.0], True),
])
def test_number_matching(written: str, sources: list[float], ok: bool) -> None:
    assert H.number_matches(written, sources) is ok


def test_check_claims_flags_only_unsupported() -> None:
    ev = H.Evidence()
    ev.add_object({"symbol": "FKBP5", "log2FC": 3.9622, "padj": 2.13e-26,
                   "gene_sets": [{"term": "HALLMARK_ANDROGEN_RESPONSE", "nes": 1.6871}]})
    text = ("FKBP5 is higher (log2FC +3.96, padj 2.1e-26) and leads HALLMARK_ANDROGEN_RESPONSE "
            "(NES 1.69). TSC22D3 rises 7.2-fold; HALLMARK_HYPOXIA is enriched.")
    c = H.check_claims(text, ev)
    assert c.unsupported_genes == ["TSC22D3"]
    assert c.unsupported_pathways == ["HALLMARK_HYPOXIA"]
    assert c.unsupported_numbers == ["7.2"]


def test_sentences_split_lines_and_sentences() -> None:
    s = H.sentences("A is up [run:a]. B is down.\n* C item\n\n1. D item")
    assert s == ["A is up [run:a].", "B is down.", "* C item", "1. D item"]


@pytest.mark.parametrize("text,status", [
    ("FKBP5 is strongly higher with dexamethasone.", "up"),
    ("VCAM1 is repressed by dexamethasone.", "down"),
    ("NR3C1 is statistically lower (padj<0.05) but below the |log2FC|>1 threshold.", "not_significant"),
    ("NR3C1 is not significantly changed.", "not_significant"),
    ("KLF15 has log2 fold change +4.39.", "up"),
    ("The weather is nice.", "not_mentioned"),
])
def test_claimed_status(text: str, status: str) -> None:
    gene = next((g for g in ["FKBP5", "VCAM1", "NR3C1", "KLF15"] if g in text), "FKBP5")
    assert claimed_status(text, gene) == status


@pytest.mark.parametrize("answer,ok", [
    ("KLF15 changes more with dexamethasone: it has a higher log2 fold change of 4.39 "
     "compared to PER1 (3.21).", True),
    ("KLF15 shows a larger change than PER1.", True),
    ("PER1 has a smaller change than KLF15.", True),
    ("PER1 changes more than KLF15.", False),
    ("KLF15 is up. PER1 is also up and more than doubled.", False),
])
def test_larger_lfc_check(answer: str, ok: bool) -> None:
    from eval.run_eval import score_question
    doc = {"eval_question": {"id": "q05", "category": "comparison", "question": "?",
                             "expect": {"genes": {"PER1": "up", "KLF15": "up"},
                                        "larger_lfc": "KLF15"}},
           "turns": [{"question": "?", "answer": answer, "events": [], "duration_s": 0,
                      "unknown_citations": []}], "runs": {}}
    s = score_question(doc, ref={})
    assert s.checks["KLF15 changes more than PER1"] is ok
