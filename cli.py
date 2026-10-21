"""Command-line entry point.

    python cli.py analyze              # QC gate, DE, enrichment (with run records)
    python cli.py lookup FKBP5         # one gene, with its provenance
    python cli.py reproduce <run_id>   # re-execute a run and verify outputs match
    python cli.py runs                 # list recent run records

    python cli.py ask "Is FKBP5 affected?"  # agent answer with citations

`eval` arrives in Phase 5.
"""

from __future__ import annotations

import argparse
import json
import sys

import pandas as pd

from agent.tools import Session
from pipeline import config, provenance


def _tag(rec: dict) -> str:
    return f"[run:{rec['run_id']}{' cache hit' if rec['cache_hit'] else ''}, {rec['duration_s']}s]"


def cmd_analyze(args: argparse.Namespace) -> int:
    s = Session()
    qc = s.run_qc()
    print(qc["summary"], _tag(qc))
    if not qc["result"]["passed"] and not args.override_qc:
        print("\nStopping: QC failed. Fix the flagged samples, or rerun with "
              "--override-qc to analyze anyway (not recommended).", file=sys.stderr)
        return 2
    de = s.run_de()
    print("\n" + de["summary"], _tag(de))
    for w in de["result"]["warnings"]:
        print(f"  note: {w}")
    t = pd.read_csv(provenance.absolute(de["outputs"]["table"]))
    top = t[t["significant"]].head(args.top)
    print(f"\nTop {len(top)} significant genes by adjusted p-value:")
    print(top[["symbol", "baseMean", "log2FC", "stat", "padj", "direction"]]
          .to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    gsea = s.run_enrichment()
    print("\n" + gsea["summary"], _tag(gsea))
    for p in gsea["result"]["top_up"] + gsea["result"]["top_down"]:
        print(f"  {p['term']:<45} NES={p['nes']:+.2f}  FDR={p['fdr']:.2g}  "
              f"leading edge: {', '.join(p['leading_edge'][:5])}")
    print("\nRun records:", ", ".join(f"runs/{r}.json" for r in s.runs))
    return 0


def cmd_lookup(args: argparse.Namespace) -> int:
    rec = Session().lookup_gene(args.symbol)
    print(rec["summary"], _tag(rec))
    if args.json:
        print(json.dumps(rec["result"], indent=2))
    return 0


def cmd_ask(args: argparse.Namespace) -> int:
    from pathlib import Path

    from agent.agent import new_agent
    agent = new_agent(Path(args.counts) if args.counts else None,
                      override_qc=args.override_qc, override_reason=args.reason)
    for q in args.question:
        ans = agent.ask(q)
        print(f"Q: {q}\n")
        print(ans.text)
        print(f"\n--- {ans.steps} model steps, {ans.duration_s}s; tools: "
              + ", ".join(f"{c['name']}({c['arguments'].get('symbol', '')})"
                          f"->{c['run_id'] or 'ERROR'}" for c in ans.tool_calls))
        print(f"--- citations: {', '.join(ans.citations) or 'none'}"
              + (f"; NOT FROM THIS SESSION: {', '.join(ans.unknown_citations)}"
                 if ans.unknown_citations else ""))
        print(f"--- transcript: {provenance.rel(ans.transcript_path)}\n")
    return 0


def cmd_reproduce(args: argparse.Namespace) -> int:
    rep = provenance.reproduce(args.run_id)
    if "error" in rep:
        print(f"REPRODUCE FAILED for {args.run_id}: {rep['error']}", file=sys.stderr)
        return 1
    print(f"Reproducing {rep['tool']} run {args.run_id} from scratch ({rep['duration_s']}s):")
    for c in rep["checks"]:
        mark = "OK  " if c["match"] else ("DIFF" if c.get("required", True) else "info")
        print(f"  {mark} {c['output']:<16} {c['method']}"
              + (f"  {c['detail']}" if c["detail"] else ""))
    for w in rep["warnings"]:
        print(f"  warning: {w}")
    print("REPRODUCED: outputs match." if rep["ok"] else "NOT REPRODUCED: outputs differ.")
    print(f"Report: runs/{args.run_id}.reproduce.json")
    return 0 if rep["ok"] else 1


def cmd_runs(args: argparse.Namespace) -> int:
    recs = sorted((json.loads(p.read_text()) for p in config.RUNS.glob("*.json")
                   if "." not in p.stem), key=lambda r: r["timestamp_utc"])
    for r in recs[-args.n:]:
        hit = "cache" if r["cache_hit"] else "fresh"
        print(f"{r['run_id']}  {r['timestamp_utc']}  {r['tool']:<15} {hit}  "
              f"{(r.get('summary') or '')[:70]}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Traceable RNA-seq interpretation agent")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("analyze", help="Run QC, differential expression, enrichment")
    p.add_argument("--top", type=int, default=20)
    p.add_argument("--override-qc", action="store_true")
    p.set_defaults(func=cmd_analyze)
    p = sub.add_parser("lookup", help="Look up one gene")
    p.add_argument("symbol")
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_lookup)
    p = sub.add_parser("ask", help="Ask the agent; several questions share one session")
    p.add_argument("question", nargs="+")
    p.add_argument("--counts", help="Alternative count matrix (e.g. a corrupted copy)")
    p.add_argument("--override-qc", action="store_true",
                   help="Let the agent interpret results even if QC fails (logged)")
    p.add_argument("--reason", default="", help="Why you are overriding QC")
    p.set_defaults(func=cmd_ask)
    p = sub.add_parser("reproduce", help="Re-execute a run and verify outputs")
    p.add_argument("run_id")
    p.set_defaults(func=cmd_reproduce)
    p = sub.add_parser("runs", help="List recent run records")
    p.add_argument("-n", type=int, default=20)
    p.set_defaults(func=cmd_runs)
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
