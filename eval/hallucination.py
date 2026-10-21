"""Deterministic claim extraction and grounding checks for agent answers.

An answer is split into a *findings* part and an optional "Context (not
computed here):" part. From each part we extract:

* gene symbols: upper-case tokens that are human gene symbols (NCBI), minus a
  stoplist of analysis words that collide with real symbols (e.g. NES is also
  the nestin gene).
* pathway names: MSigDB set names, with or without the HALLMARK_ prefix.
* numbers: after removing citations, symbols, set names, "log2", and list
  markers.

A claim is *supported* if it appears in what the agent was shown in that
session (tool results), the question, or the system prompt. Numbers match with
rounding tolerance: a value written with d decimals matches if it is within
half a unit of the d-th decimal of a source value; scientific notation matches
within 5% relative. Signs are compared on magnitude so "lower by 3.65" matches
-3.65.

Known blind spots (documented in the report): pathway names written as plain
English ("the adipogenesis pathway") and spelled-out numbers ("thirteen") are
not extracted, so they can neither pass nor fail.
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Any, Iterable

from pipeline import config, data
from pipeline.enrichment import read_gmt

CONTEXT_MARKER = re.compile(r"context\s*\(not computed here\)\s*:?", re.I)
CITATION = re.compile(r"\[run:([0-9a-zA-Z]+)\]")
SET_NAME = re.compile(r"\b(?:HALLMARK_)?[A-Z0-9]+(?:_[A-Z0-9]+)+\b")
TOKEN = re.compile(r"(?<![\w-])[A-Z][A-Z0-9-]*[A-Z0-9](?![\w-])|(?<![\w-])[A-Z](?![\w-])")
NUMBER = re.compile(r"(?<![\w.])[-+−]?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?:\s?[eE][-+−]?\d+)?(?![\w])")
SCI_TIMES = re.compile(r"(\d+(?:\.\d+)?)\s*[×x]\s*10\^?\s*([-+−]?\d+)")
LIST_MARKER = re.compile(r"^\s*(?:[-*•]|\d+[.)])\s+", re.M)

# Words that are also HGNC symbols but, in these answers, mean the analysis term.
STOPLIST = {
    "QC", "DE", "PCA", "GSEA", "NES", "FDR", "MAD", "CPM", "LFC", "ORA", "BH", "ID", "IDS",
    "RNA", "DNA", "GEO", "NCBI", "GSE", "MSIGDB", "OK", "NA", "SE", "SD", "PC", "PC1", "PC2",
    "UV", "TNF", "TNFA", "NFKB", "IL2", "STAT5", "DN", "UP", "MTORC1", "KRAS", "P53", "E2F",
    "G2M", "MYC", "WNT", "TGF", "IFN", "II", "DEX", "CELL", "SET", "FC", "MLE", "WALD", "I",
    "A", "T", "PER", "AND", "OR", "NOT", "NO", "IN", "AT", "AS", "IS", "GC", "MAX", "WAS",
}


@lru_cache(maxsize=1)
def gene_symbols() -> frozenset[str]:
    """All human gene symbols known to the pipeline's gene map."""
    gm = data.load_gene_map()
    return frozenset(gm["symbol"].dropna().astype(str))


@lru_cache(maxsize=1)
def set_names() -> frozenset[str]:
    return frozenset(read_gmt(config.GENE_SETS_FILE))


def split_context(text: str) -> tuple[str, str]:
    """(findings, context) split at the 'Context (not computed here):' marker."""
    m = CONTEXT_MARKER.search(text)
    return (text, "") if not m else (text[: m.start()], text[m.end():])


def normalize_set(name: str) -> str | None:
    """Map 'TNFA_SIGNALING_VIA_NFKB' or 'HALLMARK_TNFA...' to the GMT name, if any."""
    n = name.upper()
    if n in set_names():
        return n
    if f"HALLMARK_{n}" in set_names():
        return f"HALLMARK_{n}"
    return None


def extract_pathways(text: str) -> list[str]:
    out = []
    for m in SET_NAME.finditer(text):
        n = normalize_set(m.group(0))
        if n and n not in out:
            out.append(n)
    return out


def extract_genes(text: str) -> list[str]:
    """Gene symbols in order of first appearance (set names and citations removed)."""
    clean = SET_NAME.sub(" ", CITATION.sub(" ", text))
    syms, out = gene_symbols(), []
    for m in TOKEN.finditer(clean):
        t = m.group(0)
        if t in syms and t not in STOPLIST and t not in out:
            out.append(t)
    return out


def _to_float(s: str) -> float:
    return float(s.replace(",", "").replace("−", "-").replace(" ", ""))


def extract_numbers(text: str) -> list[str]:
    """Numeric literals, excluding citations, symbols, set names, log2/log10, list markers."""
    t = CITATION.sub(" ", text)
    t = SCI_TIMES.sub(lambda m: f"{m.group(1)}e{m.group(2)}", t)
    t = SET_NAME.sub(" ", t)
    t = re.sub(r"(?i)\b-?log(?:2|10)\b", " ", t)
    t = re.sub(r"(?i)\bpc[12]\b", " ", t)
    t = re.sub(r"\b[A-Za-z]+[0-9][A-Za-z0-9-]*\b", " ", t)  # KLF15, GSE52778, PER1
    t = re.sub(r"\b[A-Z]{2,}-\d+\b", " ", t)  # NKX2-1 style
    t = LIST_MARKER.sub(" ", t)
    return [m.group(0) for m in NUMBER.finditer(t)]


def _decimals(s: str) -> int:
    mant = re.split(r"[eE]", s.replace(",", ""))[0]
    return len(mant.split(".")[1]) if "." in mant else 0


def number_matches(written: str, sources: Iterable[float]) -> bool:
    """Does a written number match any source value within rounding tolerance?"""
    a = abs(_to_float(written))
    is_sci = bool(re.search(r"[eE]", written))
    tol = 0.5 * 10 ** (-_decimals(written)) + 1e-12
    for v in sources:
        b = abs(v)
        if is_sci or (0 < a < 1e-3):
            if b > 0 and abs(a - b) / b <= 0.05:
                return True
        elif abs(a - b) <= tol or (b != 0 and a != 0 and _decimals(written) == 0
                                   and abs(a - b) / b < 1e-9):
            return True
    return False


def _walk(obj: Any, nums: list[float], strings: list[str]) -> None:
    if isinstance(obj, dict):
        for k, v in obj.items():
            strings.append(str(k))
            _walk(v, nums, strings)
    elif isinstance(obj, list):
        for v in obj:
            _walk(v, nums, strings)
    elif isinstance(obj, bool) or obj is None:
        return
    elif isinstance(obj, (int, float)):
        if not (isinstance(obj, float) and math.isnan(obj)):
            nums.append(float(obj))
    else:
        strings.append(str(obj))


@dataclass
class Evidence:
    """What an answer may draw on: numbers, gene symbols, and set names."""

    numbers: list[float] = field(default_factory=list)
    genes: set[str] = field(default_factory=set)
    pathways: set[str] = field(default_factory=set)

    def add_object(self, obj: Any) -> None:
        nums: list[float] = []
        strings: list[str] = []
        _walk(obj, nums, strings)
        self.numbers += nums
        for s in strings:
            self.add_text(s)

    def add_text(self, text: str) -> None:
        self.numbers += [_to_float(n) for n in extract_numbers(text)]
        self.genes |= {t for t in re.findall(r"[A-Za-z0-9-]+", text) if t in gene_symbols()}
        self.pathways |= set(extract_pathways(text))


@dataclass
class ClaimCheck:
    genes: list[str]
    pathways: list[str]
    numbers: list[str]
    unsupported_genes: list[str]
    unsupported_pathways: list[str]
    unsupported_numbers: list[str]

    @property
    def n_claims(self) -> int:
        return len(self.genes) + len(self.pathways) + len(self.numbers)

    @property
    def n_unsupported(self) -> int:
        return (len(self.unsupported_genes) + len(self.unsupported_pathways)
                + len(self.unsupported_numbers))


def check_claims(text: str, evidence: Evidence) -> ClaimCheck:
    """Extract claims from text and list those not supported by the evidence."""
    genes, paths, nums = extract_genes(text), extract_pathways(text), extract_numbers(text)
    return ClaimCheck(
        genes, paths, nums,
        [g for g in genes if g not in evidence.genes],
        [p for p in paths if p not in evidence.pathways],
        [n for n in nums if not number_matches(n, evidence.numbers)],
    )


def sentences(text: str) -> list[str]:
    """Split findings into sentences; each list item counts as its own sentence."""
    out = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        m = LIST_MARKER.match(line)
        marker, rest = (m.group(0), line[m.end():]) if m else ("", line)
        parts = [s.strip() for s in re.split(r"(?<=[.!?])\s+(?=[A-Z*(])", rest) if s.strip()]
        if parts:
            parts[0] = marker + parts[0]
        out += parts
    return out


def evidence_from_session(transcript: dict[str, Any], turn: dict[str, Any]) -> Evidence:
    """Everything the agent was shown for this turn, plus question and system prompt."""
    from agent.agent import SYSTEM_PROMPT

    ev = Evidence()
    for t in transcript["turns"][: transcript["turns"].index(turn) + 1]:
        for e in t["events"]:
            if e["type"] == "tool":
                ev.add_object(e["result"])
        ev.add_text(t["question"])
    ev.add_text(SYSTEM_PROMPT)
    return ev


def evidence_from_runs(transcript: dict[str, Any], run_ids: list[str],
                       mentioned_genes: list[str], mentioned_sets: list[str]) -> Evidence:
    """Evidence contained in specific cited runs (results plus relevant table rows)."""
    import pandas as pd

    from pipeline import provenance

    ev = Evidence()
    for rid in run_ids:
        rec = transcript["runs"].get(rid)
        if rec is None:
            continue
        ev.add_object(rec["result"])
        table = rec.get("outputs", {}).get("table")
        if not table:
            continue
        path = provenance.absolute(table)
        if not path.exists():
            continue
        t = pd.read_csv(path)
        if rec["tool"] == "run_de":
            rows = t[t["symbol"].isin(mentioned_genes)]
        elif rec["tool"] == "run_enrichment":
            rows = t[t["term"].isin(mentioned_sets)]
        else:
            continue
        ev.add_object(json.loads(rows.to_json(orient="records")))
    return ev
