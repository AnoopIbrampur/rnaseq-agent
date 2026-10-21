"""Pathway enrichment tool: prerank GSEA on the DESeq2 Wald statistic.

Biology / statistics notes
--------------------------
* A **gene set** is a list of genes known to act together (here, MSigDB
  Hallmark: 50 curated processes such as "inflammatory response" or
  "adipogenesis"). Enrichment asks whether a set's genes moved together
  more than chance would predict.
* **Over-representation analysis (ORA)** takes only the significant genes,
  e.g. our 458 up-regulated, and asks with a hypergeometric test whether a set
  is overrepresented among them. It depends on an arbitrary cutoff and
  ignores the thousands of genes that moved modestly.
* **Prerank GSEA** uses *every* tested gene, ranked from most up (large
  positive Wald stat) to most down. For each set it walks down the list and
  measures whether the set's genes pile up at the top or bottom (the
  enrichment score). No cutoff is needed, and a coordinated small shift
  across many genes in a pathway can still be detected.
* **NES** (normalized enrichment score): positive = the set leans toward
  dex-induced genes, negative = toward dex-repressed genes. Significance comes
  from permutations of gene labels, and the FDR q-value corrects for testing
  50 sets.
* **Leading edge**: the genes in the set that drive the score, the ones ranked
  before the running score peaks. Those are the concrete genes to cite.
* We rank by the Wald stat rather than log2FC because it blends effect size
  with confidence: a large but noisy fold change ranks lower.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import gseapy  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

from pipeline import config  # noqa: E402


def read_gmt(path: Path) -> dict[str, list[str]]:
    """Parse a GMT file into {set name: [gene symbols]}."""
    sets = {}
    for line in path.read_text().splitlines():
        name, _url, *genes = line.split("\t")
        sets[name] = [g for g in genes if g]
    return sets


def collapse_to_symbols(table: pd.DataFrame) -> pd.DataFrame:
    """One row per symbol (config.SYMBOL_COLLAPSE_RULE).

    Prefer the row whose GeneID is current in NCBI (symbol_source == gene_info);
    among ties, keep the highest baseMean. Rows without a Wald stat are dropped.
    """
    t = table.dropna(subset=["stat", "symbol"]).copy()
    t["_current"] = t["symbol_source"] == "gene_info"
    t = t.sort_values(["symbol", "_current", "baseMean"], ascending=[True, False, False])
    return t.drop_duplicates("symbol").drop(columns="_current")


def rank_vector(table: pd.DataFrame) -> pd.Series:
    """Wald stats indexed by symbol, sorted descending with symbol as tie-breaker.

    gseapy warns that tied stats get an arbitrary order; sorting by symbol
    within ties makes that order deterministic.
    """
    t = collapse_to_symbols(table).sort_values(["stat", "symbol"], ascending=[False, True])
    return t.set_index("symbol")["stat"]


def plot_nes(res: pd.DataFrame, fdr: float, path: Path) -> None:
    """Horizontal bar chart of NES for gene sets passing the FDR threshold."""
    sig = res[res["fdr"] < fdr].sort_values("nes")
    fig, ax = plt.subplots(figsize=(7, max(2.5, 0.32 * len(sig) + 1)))
    colors = ["#E45756" if v > 0 else "#4C78A8" for v in sig["nes"]]
    ax.barh(sig["term"].str.replace("HALLMARK_", ""), sig["nes"], color=colors)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("NES (positive = higher with dex)")
    ax.set_title(f"Hallmark gene sets, prerank GSEA (FDR < {fdr})")
    ax.tick_params(axis="y", labelsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def run_enrichment(
    de_table: pd.DataFrame,
    gene_sets_path: Path = config.GENE_SETS_FILE,
    out_dir: Path = config.RESULTS / "enrichment",
    params: dict[str, Any] = config.GSEA_PARAMS,
) -> dict[str, Any]:
    """Prerank GSEA of the DE Wald stat against a local GMT file.

    Returns the top pathways in each direction (NES, FDR, leading-edge genes),
    counts at the FDR threshold, and artifact paths.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    rnk = rank_vector(de_table)
    pre = gseapy.prerank(
        rnk=rnk,
        gene_sets=str(gene_sets_path),
        min_size=params["min_size"],
        max_size=params["max_size"],
        permutation_num=params["permutation_num"],
        seed=params["seed"],
        threads=config.N_CPUS,
        no_plot=True,
        outdir=None,
        verbose=False,
    )
    res = pre.res2d.rename(
        columns={"Term": "term", "ES": "es", "NES": "nes", "NOM p-val": "pvalue",
                 "FDR q-val": "fdr", "FWER p-val": "fwer", "Tag %": "tag_pct",
                 "Gene %": "gene_pct", "Lead_genes": "leading_edge"}
    ).drop(columns=["Name"])
    for col in ("es", "nes", "pvalue", "fdr", "fwer"):
        res[col] = res[col].astype(float)
    sets = read_gmt(gene_sets_path)
    res["set_size"] = res["term"].map(lambda s: len(sets[s]))
    res["n_tested_in_set"] = res["term"].map(lambda s: int(rnk.index.isin(sets[s]).sum()))
    res = res.sort_values(["fdr", "nes"], key=lambda c: c if c.name == "fdr" else -c.abs())
    res = res.reset_index(drop=True)

    paths = {"table": out_dir / "gsea_results.csv", "ranks": out_dir / "ranks.csv",
             "plot": out_dir / "gsea_nes.png"}
    res.to_csv(paths["table"], index=False)
    rnk.rename("stat").to_csv(paths["ranks"])
    fdr = params["fdr_threshold"]
    plot_nes(res, fdr, paths["plot"])

    def top(direction: int) -> list[dict[str, Any]]:
        sel = res[(res["fdr"] < fdr) & (res["nes"] * direction > 0)]
        sel = sel.sort_values("nes", ascending=direction < 0).head(params["top_n"])
        return [
            {"term": r.term, "nes": round(r.nes, 4), "fdr": r.fdr,
             "set_size": int(r.set_size),
             "leading_edge": r.leading_edge.split(";")[:15],
             "leading_edge_size": len(r.leading_edge.split(";"))}
            for r in sel.itertuples()
        ]

    n_up = int(((res["fdr"] < fdr) & (res["nes"] > 0)).sum())
    n_down = int(((res["fdr"] < fdr) & (res["nes"] < 0)).sum())
    return {
        "summary": (
            f"Prerank GSEA on {len(rnk):,} genes ranked by Wald stat against "
            f"{len(sets)} {gene_sets_path.name} sets ({len(res)} within size limits "
            f"{params['min_size']}-{params['max_size']}): {n_up} enriched among "
            f"dex-induced genes and {n_down} among dex-repressed genes at FDR<{fdr}."
        ),
        "n_ranked_genes": int(len(rnk)),
        "n_sets_tested": int(len(res)),
        "n_sig_up": n_up,
        "n_sig_down": n_down,
        "top_up": top(1),
        "top_down": top(-1),
        "artifacts": {k: str(v) for k, v in paths.items()},
    }
