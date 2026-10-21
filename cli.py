"""Command-line entry point.

    python cli.py analyze            # QC gate, then differential expression

`ask`, `reproduce`, and `eval` arrive in later phases.
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd

from pipeline import data, de, qc


def cmd_analyze(args: argparse.Namespace) -> int:
    ds = data.load_dataset()
    qc_res = qc.run_qc(ds.counts, ds.metadata)
    print(qc_res["summary"])
    if not qc_res["passed"] and not args.override_qc:
        print("\nStopping: QC failed. Fix the flagged samples, or rerun with "
              "--override-qc to analyze anyway (not recommended).", file=sys.stderr)
        return 2
    de_res = de.run_de(ds.counts, ds.metadata, ds.gene_map)
    print("\n" + de_res["summary"])
    for w in de_res["warnings"]:
        print(f"  note: {w}")
    t = pd.read_csv(de_res["artifacts"]["table"])
    top = t[t["significant"]].head(args.top)
    print(f"\nTop {len(top)} significant genes by adjusted p-value:")
    print(top[["symbol", "baseMean", "log2FC", "stat", "padj", "direction"]]
          .to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print("\nArtifacts:")
    for path in list(qc_res["artifacts"].values()) + list(de_res["artifacts"].values()):
        print(f"  {path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Traceable RNA-seq interpretation agent")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("analyze", help="Run QC and differential expression")
    p.add_argument("--top", type=int, default=20)
    p.add_argument("--override-qc", action="store_true")
    p.set_defaults(func=cmd_analyze)
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
