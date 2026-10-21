"""Evaluation harness: run the agent on the question set and score the answers.

    python -m eval.run_eval --live          # ask the local model, save transcripts, score
    python -m eval.run_eval                 # re-score the latest saved transcripts (no LLM)
    python -m eval.run_eval --transcripts eval/transcripts/<dir>

All scoring is deterministic. Expected values come from this session's own
tool outputs (cached pipeline runs), never from hard-coded numbers.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from agent.tools import Session
from eval import hallucination as H
from pipeline import config, data, provenance

EVAL_DIR = Path(__file__).resolve().parent
TRANSCRIPT_ROOT = EVAL_DIR / "transcripts"

UP = re.compile(r"\b(higher|up\b|upregulat\w*|up-regulat\w*|increas\w*|induc\w*|elevat\w*|rais\w*)", re.I)
DOWN = re.compile(r"\b(lower|down\b|downregulat\w*|down-regulat\w*|decreas\w*|repress\w*|reduc\w*)", re.I)
NOT_SIG = re.compile(
    r"not (statistically )?significant|no significant|not counted as significant|"
    r"below the [^.]*threshold|does not meet|doesn't meet|not meet|falls short", re.I)
REFUSE_UNKNOWN = re.compile(
    r"not (a )?(known|recogni[sz]ed|valid)|unknown|not found|no (such|record)|does(n't| not) exist|"
    r"not in (this|the) (data|dataset|analysis)|don't know|do not know|cannot find|couldn't find|"
    r"not a human gene", re.I)
REFUSE_NOT_TESTED = re.compile(r"not tested|filtered|low[- ]count|too few|not enough reads|"
                               r"insufficient", re.I)
QC_BLOCK = re.compile(r"\b(qc|quality control)\b[^.]*\b(fail\w*|flag\w*|block\w*)|"
                      r"\b(fail\w*|flag\w*|block\w*)\b[^.]*\b(qc|quality control)\b", re.I)
EXPRESSION_TOOLS = {"run_de", "run_enrichment", "lookup_gene"}
WORD_NUMBERS = {w: i for i, w in enumerate(
    "zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen "
    "fifteen sixteen seventeen eighteen nineteen twenty".split())}


def cited_units(findings: str) -> list[tuple[str, list[str], bool]]:
    """(sentence, citations, inherited) for each sentence of the findings.

    List items without their own citation inherit one from the list's lead-in
    line (a cited line ending in ':') or from a citation-only line right after
    the list. `inherited` marks those, so strict coverage can also be reported.
    """
    units: list[list] = []
    block_cite: list[str] = []
    pending_items: list[list] = []
    for line in findings.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        own = H.CITATION.findall(stripped)
        is_item = bool(H.LIST_MARKER.match(line))
        if not H.CITATION.sub("", stripped).strip(" .;,") and own:  # citation-only line
            for u in pending_items:
                if not u[1]:
                    u[1], u[2] = own, True
            pending_items = []
            continue
        for sent in H.sentences(line):
            u = [sent, H.CITATION.findall(sent), False]
            if is_item:
                if not u[1] and block_cite:
                    u[1], u[2] = list(block_cite), True
                pending_items.append(u)
            units.append(u)
        if not is_item:
            block_cite = own if stripped.rstrip("*").endswith(":") else []
            pending_items = [] if not stripped.endswith(":") else pending_items
    return [tuple(u) for u in units]


# --- reference values from the pipeline -------------------------------------------


def reference() -> dict[str, Any]:
    """Pipeline outputs on the real data (cached), used to resolve expectations."""
    s = Session()
    de, gsea = s.run_de(), s.run_enrichment()
    table = pd.read_csv(provenance.absolute(de["outputs"]["table"]))
    return {"run_de": de["result"], "run_enrichment": gsea["result"], "de_table": table,
            "de_run_id": de["run_id"], "gsea_run_id": gsea["run_id"]}


def resolve(ref: dict[str, Any], path: str) -> Any:
    tool, key = path.split(".", 1)
    return ref[tool][key]


def verify_expected(ref: dict[str, Any]) -> list[dict[str, Any]]:
    """Check expected.yaml genes against the actual DE results."""
    exp = yaml.safe_load((EVAL_DIR / "expected.yaml").read_text())["genes"]
    t = ref["de_table"].set_index("symbol")
    rows = []
    for g, direction in exp.items():
        if g not in t.index:
            rows.append({"gene": g, "expected": direction, "observed": "not tested",
                         "log2FC": None, "padj": None, "ok": False})
            continue
        r = t.loc[g] if not isinstance(t.loc[g], pd.DataFrame) else t.loc[g].iloc[0]
        rows.append({"gene": g, "expected": direction, "observed": r["direction"],
                     "log2FC": round(float(r["log2FC"]), 3), "padj": float(r["padj"]),
                     "ok": r["direction"] == direction})
    return rows


def authors_concordance(ref: dict[str, Any]) -> dict[str, Any]:
    """Compare our DE with the study authors' Cuffdiff results (independent pipeline)."""
    a = pd.read_csv(config.AUTHORS_DE_FILE, sep="\t")
    a = a[a["status"] == "OK"].drop_duplicates("gene").set_index("gene")
    # Cuffdiff's sample_1=Dex, sample_2=Untreated, log2(fold_change) = log2(value_2/value_1)
    lfc = pd.to_numeric(a["log2(fold_change)"], errors="coerce").replace([float("inf"), float("-inf")], pd.NA)
    sign_check = a[["sample_1", "sample_2"]].drop_duplicates().to_dict(orient="records")
    a_dex_up = -lfc if sign_check[0]["sample_1"] == "Dex" else lfc
    ours = ref["de_table"].drop_duplicates("symbol").set_index("symbol")
    common = ours.index.intersection(a.index)
    both = [g for g in common if ours.loc[g, "significant"] and a.loc[g, "significant"] == "yes"]
    agree = sum((ours.loc[g, "log2FC"] > 0) == (a_dex_up[g] > 0) for g in both if pd.notna(a_dex_up[g]))
    ours_sig = set(ours.index[ours["significant"]])
    theirs_sig = set(a.index[a["significant"] == "yes"])
    corr = pd.concat([ours.loc[common, "log2FC_mle"], a_dex_up[common]], axis=1).dropna().astype(float)
    return {
        "cuffdiff_columns": sign_check[0],
        "genes_in_common": int(len(common)),
        "sig_ours": len(ours_sig), "sig_authors": len(theirs_sig),
        "sig_both": len(both), "direction_agree_in_both": int(agree),
        "jaccard_sig": round(len(ours_sig & theirs_sig) / max(1, len(ours_sig | theirs_sig)), 3),
        "spearman_lfc": round(float(corr.corr(method="spearman").iloc[0, 1]), 3),
    }


# --- live run ------------------------------------------------------------------


def run_live(questions: list[dict[str, Any]], out_dir: Path) -> None:
    """Ask each question in a fresh agent session and save one transcript per question."""
    from agent.agent import new_agent
    out_dir.mkdir(parents=True, exist_ok=True)
    corrupted = data.write_corrupted_counts("scrambled_sample")
    for q in questions:
        agent = new_agent(corrupted if q.get("data") == "corrupted" else None)
        t0 = datetime.now()
        try:
            ans = agent.ask(q["question"], save=False)
            print(f"{q['id']} ({ans.duration_s}s): {ans.text[:90]!r}", flush=True)
        except RuntimeError as e:
            print(f"{q['id']} FAILED: {e}", file=sys.stderr)
            raise
        path = agent.save_transcript(out_dir / f"{q['id']}.json")
        doc = json.loads(path.read_text())
        doc["eval_question"] = q
        doc["started"] = t0.isoformat(timespec="seconds")
        path.write_text(json.dumps(doc, indent=2, default=str))


# --- scoring -----------------------------------------------------------------------


def claimed_status(findings: str, gene: str) -> str:
    """What the answer claims about one gene: up, down, not_significant, or unclear."""
    sents = [s for s in H.sentences(findings) if re.search(rf"\b{re.escape(gene)}\b", s)]
    if not sents:
        return "not_mentioned"
    text = " ".join(sents)
    if NOT_SIG.search(text):
        return "not_significant"
    up, down = bool(UP.search(text)), bool(DOWN.search(text))
    if up and not down:
        return "up"
    if down and not up:
        return "down"
    m = re.search(rf"\b{re.escape(gene)}\b[^.]*?log2\s*(?:fold change|fc)?[^-+\d]*([-+−]?\d+(?:\.\d+)?)", text, re.I)
    if m:
        return "up" if H._to_float(m.group(1)) > 0 else "down"
    return "unclear"


@dataclass
class QuestionScore:
    id: str
    category: str
    question: str
    answer: str
    duration_s: float
    tools: list[str]
    checks: dict[str, bool] = field(default_factory=dict)
    gene_claims: dict[str, dict[str, str]] = field(default_factory=dict)
    n_claims: int = 0
    unsupported: dict[str, list[str]] = field(default_factory=dict)
    context_violations: list[str] = field(default_factory=list)
    factual_sentences: int = 0
    cited_sentences: int = 0
    strictly_cited_sentences: int = 0
    supported_cited_sentences: int = 0
    unsupported_citation_details: list[str] = field(default_factory=list)
    unknown_citations: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return all(self.checks.values()) and not any(self.unsupported.values()) \
            and not self.unknown_citations


def score_question(doc: dict[str, Any], ref: dict[str, Any]) -> QuestionScore:
    q = doc["eval_question"]
    turn = doc["turns"][-1]
    answer = turn["answer"]
    findings, context = H.split_context(answer)
    tools = [f"{e['name']}{'!' if 'error' in e['result'] else ''}"
             for e in turn["events"] if e["type"] == "tool"]
    s = QuestionScore(q["id"], q["category"], q["question"], answer, turn["duration_s"], tools,
                      unknown_citations=turn["unknown_citations"])
    exp = q.get("expect", {})

    # Grounding: every gene/pathway/number vs everything the agent was shown.
    ev = H.evidence_from_session(doc, turn)
    chk = H.check_claims(findings, ev)
    s.n_claims = chk.n_claims
    s.unsupported = {"genes": chk.unsupported_genes, "pathways": chk.unsupported_pathways,
                     "numbers": chk.unsupported_numbers}
    s.context_violations = H.extract_numbers(context) + H.CITATION.findall(context)

    # Citation coverage and support (claims in a cited sentence must be in the cited runs).
    for sent, cited, inherited in cited_units(findings):
        c = H.check_claims(sent, ev)
        if c.n_claims == 0:
            continue
        s.factual_sentences += 1
        if not cited:
            continue
        s.cited_sentences += 1
        s.strictly_cited_sentences += not inherited
        run_ev = H.evidence_from_runs(doc, cited, c.genes, c.pathways)
        rc = H.check_claims(sent, run_ev)
        if rc.n_unsupported == 0:
            s.supported_cited_sentences += 1
        else:
            s.unsupported_citation_details.append(
                f"{sent[:120]!r} cites {cited} but those runs lack "
                f"{rc.unsupported_genes + rc.unsupported_pathways + rc.unsupported_numbers}")

    # Expectations.
    for gene, want in exp.get("genes", {}).items():
        got = claimed_status(findings, gene)
        s.gene_claims[gene] = {"expected": want, "claimed": got}
        s.checks[f"gene {gene} {want}"] = got == want
    for path in exp.get("numbers", []):
        val = resolve(ref, path)
        written = [H._to_float(n) for n in H.extract_numbers(findings)]
        written += [v for w, v in WORD_NUMBERS.items() if re.search(rf"\b{w}\b", findings, re.I)]
        ok = H.number_matches(str(val), written) if written else False
        if not ok and val == 0:
            ok = bool(re.search(r"\b(no|none|zero|not any)\b", findings, re.I))
        s.checks[f"number {path}={val}"] = ok
    for p in exp.get("pathways", []):
        name = p.replace("HALLMARK_", "").replace("_", " ")
        s.checks[f"pathway {p}"] = p in H.extract_pathways(findings) or name.lower() in findings.lower()
    if "order" in exp:
        n = exp.get("order_n", 5)
        want = [g["symbol"] for g in resolve(ref, exp["order"])[:n]]
        got = [g for g in H.extract_genes(findings) if g in want]
        s.checks[f"ranking top {n}"] = got == want
        if got != want:
            s.notes.append(f"expected order {want}, answer order {got}")
    if "larger_lfc" in exp:
        big = exp["larger_lfc"]
        small = next(g for g in exp["genes"] if g != big)
        cmp = r"(?:than|compared (?:to|with)|vs\.?|versus|relative to)"
        ok = any(
            re.search(rf"\b{big}\b.*\b(more|larger|greater|stronger|bigger|higher)\b.*"
                      rf"\b{cmp}\b.*\b{small}\b", sent, re.I)
            or re.search(rf"\b{small}\b.*\b(less|smaller|weaker|lower)\b.*\b{cmp}\b.*\b{big}\b",
                         sent, re.I)
            for sent in H.sentences(findings))
        s.checks[f"{big} changes more than {small}"] = ok
    refusal = exp.get("refusal")
    if refusal == "unknown_symbol":
        s.checks["refuses unknown gene"] = bool(REFUSE_UNKNOWN.search(findings)) and \
            not any(chk.unsupported_numbers)
    elif refusal == "not_tested":
        s.checks["says not tested"] = bool(REFUSE_NOT_TESTED.search(findings))
    elif refusal == "qc_blocked":
        ran = [e for e in turn["events"] if e["type"] == "tool" and e["name"] in EXPRESSION_TOOLS
               and "error" not in e["result"]]
        s.checks["no expression tool succeeded"] = not ran
        s.checks["answer explains QC failure"] = bool(QC_BLOCK.search(findings))
    if exp.get("no_genes"):
        s.checks["names no genes"] = not H.extract_genes(findings)
    return s


def summarize(scores: list[QuestionScore], verified: list[dict[str, Any]]) -> dict[str, Any]:
    verified_genes = {r["gene"] for r in verified if r["ok"]}
    agree = [(g, c) for s in scores for g, c in s.gene_claims.items() if g in verified_genes]
    refusal = [v for s in scores for k, v in s.checks.items()
               if s.category in ("trick", "qc_gate")]
    claims = sum(s.n_claims for s in scores)
    unsup = sum(len(v) for s in scores for v in s.unsupported.values())
    return {
        "questions": len(scores),
        "passed": sum(s.passed for s in scores),
        "agreement": (sum(c["claimed"] == c["expected"] for _, c in agree), len(agree)),
        "answers_with_hallucination": sum(any(s.unsupported.values()) for s in scores),
        "claims_total": claims,
        "claims_unsupported": unsup,
        "citation_coverage": (sum(s.cited_sentences for s in scores),
                              sum(s.factual_sentences for s in scores)),
        "citation_coverage_strict": (sum(s.strictly_cited_sentences for s in scores),
                                     sum(s.factual_sentences for s in scores)),
        "citation_support": (sum(s.supported_cited_sentences for s in scores),
                             sum(s.cited_sentences for s in scores)),
        "fabricated_citations": sum(len(s.unknown_citations) for s in scores),
        "refusal_checks": (sum(refusal), len(refusal)),
        "context_violations": sum(len(s.context_violations) for s in scores),
        "mean_latency_s": round(sum(s.duration_s for s in scores) / max(1, len(scores)), 1),
    }


def pct(pair: tuple[int, int]) -> str:
    n, d = pair
    return f"{n}/{d} ({n / d:.0%})" if d else "n/a"


def write_report(scores: list[QuestionScore], summ: dict[str, Any], verified: list[dict[str, Any]],
                 conc: dict[str, Any], tdir: Path, meta: dict[str, Any]) -> str:
    L = ["# Evaluation report", "",
         f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')} from transcripts in "
         f"`{provenance.rel(tdir)}`. Model: `{meta.get('model')}` (local Ollama, think={meta.get('think')}, "
         f"temperature 0). System prompt sha256 `{(meta.get('prompt_sha') or '')[:12]}`.", "",
         "All checks are deterministic code (`eval/hallucination.py`, `eval/run_eval.py`); no LLM "
         "grades another LLM. Expected values are read from the pipeline's own run outputs.", "",
         "## Headline numbers", "",
         "| Metric | Result | What it means |", "|---|---|---|",
         f"| Questions fully passed | {summ['passed']}/{summ['questions']} | every expectation met, no unsupported claim, no fabricated citation |",
         f"| Agreement with expected genes | {pct(summ['agreement'])} | literature genes reported as significant in the right direction |",
         f"| Answers with ≥1 unsupported claim | {summ['answers_with_hallucination']}/{summ['questions']} | a gene, pathway, or number not in that session's tool outputs |",
         f"| Unsupported claims / all claims | {summ['claims_unsupported']}/{summ['claims_total']} | claim-level hallucination rate |",
         f"| Citation coverage | {pct(summ['citation_coverage'])} | factual sentences carrying a `[run:…]` citation (list items may inherit the list's citation) |",
         f"| Citation coverage, strict | {pct(summ['citation_coverage_strict'])} | same, but every sentence and list item needs its own citation |",
         f"| Citation support | {pct(summ['citation_support'])} | cited sentences whose claims are in the cited runs |",
         f"| Fabricated run ids | {summ['fabricated_citations']} | citations to runs that don't exist in the session |",
         f"| Refusal / QC-gate checks | {pct(summ['refusal_checks'])} | unknown gene, untested gene, failed QC |",
         f"| Numbers or citations inside 'Context' | {summ['context_violations']} | rule: context paragraphs carry no numbers |",
         f"| Mean answer latency | {summ['mean_latency_s']} s | local 9B model on a laptop |", "",
         "## Per question", "",
         "| id | category | pass | checks | unsupported claims | cited/factual sentences | tools | s |",
         "|---|---|---|---|---|---|---|---|"]
    for s in scores:
        unsup = ", ".join(v for vs in s.unsupported.values() for v in vs) or "–"
        checks = "; ".join(f"{'✅' if v else '❌'} {k}" for k, v in s.checks.items()) or "–"
        L.append(f"| {s.id} | {s.category} | {'✅' if s.passed else '❌'} | {checks} | {unsup} | "
                 f"{s.cited_sentences}/{s.factual_sentences} | {', '.join(s.tools) or '–'} | {s.duration_s} |")
    L += ["", "## Failures and issues (every one)", ""]
    any_issue = False
    for s in scores:
        issues = [f"check failed: {k}" for k, v in s.checks.items() if not v]
        issues += [f"unsupported {kind}: {', '.join(v)}" for kind, v in s.unsupported.items() if v]
        issues += [f"fabricated citation: {c}" for c in s.unknown_citations]
        issues += [f"uncited factual sentences: {s.factual_sentences - s.cited_sentences}"] \
            if s.cited_sentences < s.factual_sentences else []
        issues += [f"citation does not support claim: {d}" for d in s.unsupported_citation_details]
        issues += [f"context paragraph contains: {', '.join(s.context_violations)}"] if s.context_violations else []
        issues += s.notes
        if not issues:
            continue
        any_issue = True
        L += [f"### {s.id}: {s.question}", ""] + [f"- {i}" for i in issues] + \
             ["", "<details><summary>Answer</summary>", "", "```", s.answer, "```", "", "</details>", ""]
    if not any_issue:
        L.append("None.")
    L += ["", "## Pipeline sanity checks (no LLM involved)", "",
          "### Expected genes vs this pipeline's DE results", "",
          "Sanity check, not ground truth: see `eval/expected.yaml` for where the list comes from.", "",
          "| gene | expected | observed | log2FC (shrunken) | padj | ok |", "|---|---|---|---|---|---|"]
    for r in verified:
        padj = f"{r['padj']:.2g}" if r["padj"] is not None else "–"
        L.append(f"| {r['gene']} | {r['expected']} | {r['observed']} | {r['log2FC']} | {padj} | {'✅' if r['ok'] else '❌'} |")
    L += ["", "### Agreement with the study authors' own analysis", "",
          "The authors' Cuffdiff dex-vs-untreated table from GEO "
          "(`GSE52778_Dex_vs_Untreated_gene_exp.diff.gz`) is an independent pipeline (different aligner, "
          "annotation, statistical model, and no paired design), so perfect overlap is not expected.", "",
          f"- Genes matched by symbol: {conc['genes_in_common']:,}",
          f"- Significant here: {conc['sig_ours']:,} (padj<0.05, |log2FC|>1); significant in Cuffdiff: {conc['sig_authors']:,} (its own q<0.05 call)",
          f"- Significant in both: {conc['sig_both']:,}; same direction: {conc['direction_agree_in_both']:,}",
          f"- Jaccard overlap of significant sets: {conc['jaccard_sig']}",
          f"- Spearman correlation of log2 fold changes (unshrunken here): {conc['spearman_lfc']}", "",
          "## Method notes and known blind spots", "",
          "- **Hallucination** = a gene symbol, MSigDB set name, or number in the findings part of an answer "
          "that does not appear in anything the agent was shown in that session (tool results, the question, "
          "the system prompt). Numbers match with rounding tolerance (half a unit in the last written decimal; "
          "5% relative for scientific notation).",
          "- **Not detected:** pathway names written in plain English, spelled-out numbers (\"thirteen\"), "
          "and wrong *relationships* between correct facts (e.g. right number, wrong gene). The ranking check "
          "covers one such case explicitly.",
          "- **Direction claims** are read with keyword rules (higher/up/induced vs lower/down/repressed, plus "
          "'not significant' phrases). Ambiguous sentences are scored `unclear`, which counts as wrong.",
          "- **Citation support** re-checks each cited sentence against only the runs it cites (their results "
          "plus the relevant rows of their output tables).",
          "- One run of 15 questions at temperature 0 is a small sample; rerunning can change answers. "
          "See `eval/CHANGELOG.md` for any prompt changes made after looking at results.", ""]
    return "\n".join(L)


def latest_transcript_dir() -> Path:
    dirs = sorted(p for p in TRANSCRIPT_ROOT.iterdir() if p.is_dir() and p.name.startswith("run_"))
    if not dirs:
        raise FileNotFoundError(
            f"No saved eval transcripts in {TRANSCRIPT_ROOT}. Run with --live first "
            "(needs Ollama running with the model pulled).")
    return dirs[-1]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--live", action="store_true", help="ask the local model (else re-score saved)")
    ap.add_argument("--transcripts", type=Path, help="score this transcript directory")
    ap.add_argument("--only", nargs="*", help="question ids to run")
    args = ap.parse_args(argv)

    questions = yaml.safe_load((EVAL_DIR / "questions.yaml").read_text())["questions"]
    if args.only:
        questions = [q for q in questions if q["id"] in args.only]
    if args.live:
        tdir = TRANSCRIPT_ROOT / f"run_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"
        run_live(questions, tdir)
    else:
        tdir = args.transcripts or latest_transcript_dir()
        print(f"Scoring saved transcripts in {provenance.rel(tdir)} (no model calls).")

    ref = reference()
    docs = [json.loads(p.read_text()) for p in sorted(tdir.glob("q*.json"))]
    scores = [score_question(d, ref) for d in docs]
    verified = verify_expected(ref)
    conc = authors_concordance(ref)
    summ = summarize(scores, verified)
    meta = {"model": docs[0]["backend"]["model"], "think": docs[0]["backend"]["think"],
            "prompt_sha": docs[0]["system_prompt_sha256"]} if docs else {}
    report = write_report(scores, summ, verified, conc, tdir, meta)
    (EVAL_DIR / "report.md").write_text(report)
    (tdir / "scores.json").write_text(json.dumps(
        {"summary": summ, "scores": [asdict(s) | {"passed": s.passed} for s in scores],
         "expected_verified": verified, "authors_concordance": conc}, indent=2, default=str))
    print(f"\nPassed {summ['passed']}/{summ['questions']}; agreement {pct(summ['agreement'])}; "
          f"answers with unsupported claims {summ['answers_with_hallucination']}; "
          f"citation coverage {pct(summ['citation_coverage'])}; refusal {pct(summ['refusal_checks'])}.")
    print(f"Report: {provenance.rel(EVAL_DIR / 'report.md')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
