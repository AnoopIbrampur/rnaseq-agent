"""Differential expression tool: PyDESeq2, paired design, dex vs untreated.

Biology / statistics notes
--------------------------
* **Design `~ cell + dex`.** Each donor cell line gave one untreated and one
  treated sample. Donors differ a lot at baseline (they separate on PC2 in
  QC). Including `cell` as a *blocking factor* lets the model estimate the dex
  effect *within* each donor, so donor differences don't count as noise. This
  is the RNA-seq version of a paired t-test.
* **Normalization (size factors).** DESeq2 estimates one scaling factor per
  sample so that a typical gene has the same expression across samples. It
  corrects for sequencing depth without changing the raw counts it models.
* **Dispersion.** Read counts are noisier than a Poisson model predicts, and
  how much noisier differs by gene. DESeq2 estimates this extra variance
  (dispersion) per gene, borrowing strength across genes because 8 samples is
  too few to estimate it gene by gene.
* **log2 fold change (LFC).** log2(treated / untreated). +1 means doubled with
  dex, -1 halved. **Shrinkage** pulls noisy LFCs of low-count genes toward 0
  so they don't look dramatic by chance; the Wald stat and p-value are not
  changed by it.
* **Wald test + Benjamini-Hochberg.** For each gene, the Wald stat is
  LFC / standard error. With ~16k genes tested, ~800 would hit p<0.05 by luck;
  BH-adjusted p-values (padj) control the false discovery rate instead.
* `padj` is NaN when DESeq2 drops a gene by *independent filtering* (too
  little signal to ever be significant) or the gene has a *Cook's distance*
  outlier (one sample dominates its counts).
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from pydeseq2.dds import DeseqDataSet  # noqa: E402
from pydeseq2.default_inference import DefaultInference  # noqa: E402
from pydeseq2.ds import DeseqStats  # noqa: E402

from pipeline import config  # noqa: E402
from pipeline.qc import low_count_mask  # noqa: E402


def _design_metadata(metadata: pd.DataFrame) -> pd.DataFrame:
    """Metadata with factor levels fixed so 'untrt' is the reference.

    PyDESeq2 (formulaic) uses the first category as the reference level. With
    plain strings it would sort alphabetically, making 'trt' the reference and
    flipping the sign of the shrunken LFC. Ordered categories prevent that.
    """
    meta = metadata[["cell", "dex"]].copy()
    meta["cell"] = pd.Categorical(meta["cell"].astype(str))
    meta["dex"] = pd.Categorical(meta["dex"].astype(str), categories=["untrt", "trt"])
    return meta


def fit_deseq2(
    counts: pd.DataFrame, metadata: pd.DataFrame, params: dict[str, Any]
) -> tuple[pd.DataFrame, list[str]]:
    """Fit DESeq2 and return (results indexed by GeneID, warnings raised).

    `counts` is genes x samples, already low-count filtered.
    """
    np.random.seed(config.SEED)
    meta = _design_metadata(metadata)
    inference = DefaultInference(n_cpus=config.N_CPUS)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        dds = DeseqDataSet(
            counts=counts.T.rename(columns=str),
            metadata=meta,
            design=params["design"],
            refit_cooks=True,
            inference=inference,
            quiet=True,
        )
        dds.deseq2()
        coeffs = list(dds.varm["LFC"].columns)
        if params["shrink_coeff"] not in coeffs:
            raise RuntimeError(
                f"Expected coefficient {params['shrink_coeff']!r} not in model "
                f"coefficients {coeffs}. Check the design and reference levels."
            )
        stats = DeseqStats(
            dds, contrast=params["contrast"], alpha=params["alpha"],
            inference=inference, quiet=True,
        )
        stats.summary()
        mle = stats.results_df[["log2FoldChange", "lfcSE"]].copy()
        stats.lfc_shrink(coeff=params["shrink_coeff"])
    res = stats.results_df.copy()
    res["log2FoldChange_mle"] = mle["log2FoldChange"]
    res["lfcSE_mle"] = mle["lfcSE"]
    res.index = res.index.astype(int)
    msgs = sorted({f"{w.category.__name__}: {w.message}" for w in caught})
    return res, msgs


def annotate(res: pd.DataFrame, gene_map: pd.DataFrame, params: dict[str, Any]) -> pd.DataFrame:
    """Add symbols, significance, direction, and the reason for any missing padj."""
    gm = gene_map.set_index("gene_id")
    out = pd.DataFrame(
        {
            "gene_id": res.index,
            "symbol": gm.reindex(res.index)["symbol"].to_numpy(),
            "symbol_source": gm.reindex(res.index)["symbol_source"].to_numpy(),
            "baseMean": res["baseMean"].to_numpy(),
            "log2FC": res["log2FoldChange"].to_numpy(),
            "lfcSE": res["lfcSE"].to_numpy(),
            "log2FC_mle": res["log2FoldChange_mle"].to_numpy(),
            "stat": res["stat"].to_numpy(),
            "pvalue": res["pvalue"].to_numpy(),
            "padj": res["padj"].to_numpy(),
        }
    )
    sig = (out["padj"] < params["alpha"]) & (out["log2FC"].abs() > params["lfc_threshold"])
    out["significant"] = sig
    out["direction"] = np.where(sig, np.where(out["log2FC"] > 0, "up", "down"), "ns")
    out["note"] = np.select(
        [out["pvalue"].isna(), out["padj"].isna()],
        ["cooks_outlier", "independent_filtering"],
        default="",
    )
    return out.sort_values(["padj", "pvalue"], na_position="last").reset_index(drop=True)


def plot_volcano(table: pd.DataFrame, params: dict[str, Any], path: Path, label_top: int = 12) -> None:
    """Volcano plot: shrunken log2FC vs -log10(padj), significant genes colored."""
    t = table.dropna(subset=["padj"])
    y = -np.log10(t["padj"].clip(lower=1e-300))
    colors = t["direction"].map({"up": "#E45756", "down": "#4C78A8", "ns": "#BBBBBB"})
    fig, ax = plt.subplots(figsize=(7, 5.5))
    ax.scatter(t["log2FC"], y, c=colors, s=6, alpha=0.7, linewidths=0)
    ax.axhline(-np.log10(params["alpha"]), color="gray", ls="--", lw=0.8)
    for x in (-params["lfc_threshold"], params["lfc_threshold"]):
        ax.axvline(x, color="gray", ls="--", lw=0.8)
    for _, r in t[t["significant"]].head(label_top).iterrows():
        ax.annotate(r["symbol"], (r["log2FC"], -np.log10(max(r["padj"], 1e-300))),
                    fontsize=7, xytext=(3, 3), textcoords="offset points")
    n_up, n_down = (t["direction"] == "up").sum(), (t["direction"] == "down").sum()
    ax.set_xlabel("log2 fold change (dex vs untreated, shrunken)")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title(f"Dexamethasone response: {n_up} up, {n_down} down "
                 f"(padj<{params['alpha']}, |LFC|>{params['lfc_threshold']:g})")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def run_de(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    gene_map: pd.DataFrame,
    out_dir: Path = config.RESULTS / "de",
    params: dict[str, Any] = config.DE_PARAMS,
    low_count: dict[str, int] = config.LOW_COUNT_FILTER,
) -> dict[str, Any]:
    """Run differential expression for dex vs untreated, controlling for cell line.

    Uses raw counts. Returns summary numbers, the top genes, the list of genes
    removed by the low-count filter, any modelling warnings, and artifact paths.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    keep = low_count_mask(counts, **low_count)
    res, msgs = fit_deseq2(counts.loc[keep], metadata, params)
    table = annotate(res, gene_map, params)

    paths = {
        "table": out_dir / "de_results.csv",
        "filtered_genes": out_dir / "low_count_filtered.csv",
        "volcano": out_dir / "volcano.png",
    }
    table.to_csv(paths["table"], index=False)
    gm = gene_map.set_index("gene_id")
    filtered = counts.index[~keep]
    pd.DataFrame(
        {"gene_id": filtered, "symbol": gm.reindex(filtered)["symbol"].to_numpy(),
         "max_count": counts.loc[~keep].max(axis=1).to_numpy()}
    ).to_csv(paths["filtered_genes"], index=False)
    plot_volcano(table, params, paths["volcano"])

    n_up = int((table["direction"] == "up").sum())
    n_down = int((table["direction"] == "down").sum())
    tested = int(table["padj"].notna().sum())
    cols = ["gene_id", "symbol", "baseMean", "log2FC", "stat", "padj"]
    return {
        "summary": (
            f"DE (dex vs untreated, design {params['design']}): {len(table):,} genes "
            f"modelled after low-count filter, {tested:,} with an adjusted p-value. "
            f"{n_up} up and {n_down} down at padj<{params['alpha']} and "
            f"|log2FC|>{params['lfc_threshold']:g}."
        ),
        "n_genes_modelled": int(len(table)),
        "n_genes_with_padj": tested,
        "n_low_count_filtered": int((~keep).sum()),
        "n_sig_up": n_up,
        "n_sig_down": n_down,
        "n_padj_lt_alpha_any_lfc": int((table["padj"] < params["alpha"]).sum()),
        "top_up": table[table["direction"] == "up"].head(20)[cols].to_dict(orient="records"),
        "top_down": table[table["direction"] == "down"].head(20)[cols].to_dict(orient="records"),
        "warnings": msgs,
        "artifacts": {k: str(v) for k, v in paths.items()},
    }
