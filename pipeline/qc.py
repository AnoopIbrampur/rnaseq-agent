"""Quality control tool with a gate that can block interpretation.

Biology notes
-------------
* **Library size** is a sample's total read count. It varies with how much
  was sequenced, not with biology, so a very small library is a technical
  warning sign rather than a finding.
* **log-CPM** (log2 counts per million) rescales each sample to the same
  total, then log-transforms so a few huge genes don't dominate. It's fine
  for *looking* at samples (PCA, correlations) but not for the DE test.
* **PCA** squashes thousands of genes into a few axes of largest variation.
  If treated and untreated samples land in different places on PC1, the
  treatment is the strongest signal in the data, which is what we hope to see.
  Donor cell lines often separate on a later PC; that's the donor effect the
  paired design (`~cell + dex`) controls for.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from pipeline import config  # noqa: E402


def low_count_mask(counts: pd.DataFrame, min_count: int, min_samples: int) -> pd.Series:
    """True for genes with >= min_count reads in at least min_samples samples."""
    return (counts >= min_count).sum(axis=1) >= min_samples


def log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """log2(counts per million + 1), per sample."""
    return np.log2(counts / counts.sum() * 1e6 + 1)


def pca(values: pd.DataFrame, top_genes: int) -> tuple[pd.DataFrame, np.ndarray]:
    """PCA of samples using the `top_genes` most variable genes.

    `values` is genes x samples. Returns (scores: samples x PCs, explained
    variance ratio per PC). Uses SVD, so it is deterministic up to sign; signs
    are fixed so each PC's largest-magnitude loading is positive.
    """
    top = values.loc[values.var(axis=1).sort_values(ascending=False).index[:top_genes]]
    x = (top.T - top.T.mean()).to_numpy()
    u, s, vt = np.linalg.svd(x, full_matrices=False)
    signs = np.sign(vt[np.arange(len(vt)), np.abs(vt).argmax(axis=1)])
    u, vt = u * signs, vt * signs[:, None]
    scores = pd.DataFrame(
        u * s, index=values.columns, columns=[f"PC{i + 1}" for i in range(len(s))]
    )
    return scores, s**2 / (s**2).sum()


def eta_squared(scores: pd.Series, groups: pd.Series) -> float:
    """Fraction of variance in `scores` explained by `groups` (one-way ANOVA R^2)."""
    total = ((scores - scores.mean()) ** 2).sum()
    if total == 0:
        return 0.0
    between = sum(
        len(g) * (g.mean() - scores.mean()) ** 2
        for _, g in scores.groupby(groups.to_numpy(), observed=True)
    )
    return float(between / total)


def library_size_flags(lib_sizes: pd.Series, min_fraction: float) -> list[dict[str, Any]]:
    """Flag samples whose library size is below min_fraction of the median."""
    threshold = min_fraction * lib_sizes.median()
    return [
        {
            "sample": s,
            "rule": "min_library_fraction",
            "value": int(n),
            "threshold": float(threshold),
            "message": (
                f"{s} has {n / 1e6:.1f}M reads, below {min_fraction:.0%} of the "
                f"median ({lib_sizes.median() / 1e6:.1f}M)."
            ),
        }
        for s, n in lib_sizes.items()
        if n < threshold
    ]


def pca_outlier_flags(
    scores: pd.DataFrame, groups: pd.Series, n_mad: float
) -> tuple[list[dict[str, Any]], pd.Series]:
    """Flag samples far from their group's centroid on PC1/PC2.

    Returns (flags, distance per sample). The cutoff is
    median(distance) + n_mad * 1.4826 * MAD(distance), computed over all samples.
    """
    pc = scores[["PC1", "PC2"]]
    centroids = pc.groupby(groups.to_numpy(), observed=True).transform("mean")
    dist = np.sqrt(((pc - centroids) ** 2).sum(axis=1))
    med = dist.median()
    robust_sd = 1.4826 * (dist - med).abs().median()
    cutoff = med + n_mad * robust_sd
    dist.attrs["cutoff"] = float(cutoff)
    flags = [
        {
            "sample": s,
            "rule": "pca_outlier_mad",
            "value": float(d),
            "threshold": float(cutoff),
            "message": (
                f"{s} is {d:.1f} PCA units from its group centroid; cutoff is "
                f"{cutoff:.1f} (median {med:.1f} + {n_mad:g} robust SD)."
            ),
        }
        for s, d in dist.items()
        if robust_sd > 0 and d > cutoff
    ]
    return flags, dist


def separation_report(
    scores: pd.DataFrame, explained: np.ndarray, metadata: pd.DataFrame, r2: float
) -> dict[str, Any]:
    """Which factors (dex, cell) explain PC1 and PC2."""
    report: dict[str, Any] = {}
    for i, pc in enumerate(["PC1", "PC2"]):
        report[pc] = {
            "variance_explained": float(explained[i]),
            "r2_dex": eta_squared(scores[pc], metadata["dex"]),
            "r2_cell": eta_squared(scores[pc], metadata["cell"]),
        }
    by = {
        f: [pc for pc in ("PC1", "PC2") if report[pc][f"r2_{f}"] >= r2] for f in ("dex", "cell")
    }
    if by["dex"] and by["cell"]:
        verdict = "both"
    elif by["dex"]:
        verdict = "treatment"
    elif by["cell"]:
        verdict = "cell line"
    else:
        verdict = "neither"
    report["separates_by"] = verdict
    report["dex_pcs"], report["cell_pcs"] = by["dex"], by["cell"]
    return report


def _plot_library_sizes(lib: pd.Series, threshold: float, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.bar(lib.index, lib / 1e6, color="#4C78A8")
    ax.axhline(threshold / 1e6, color="#E45756", ls="--", label="QC threshold")
    ax.set_ylabel("Reads (millions)")
    ax.set_title("Library size per sample")
    ax.tick_params(axis="x", rotation=45)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def _plot_pca(
    scores: pd.DataFrame, explained: np.ndarray, metadata: pd.DataFrame,
    flagged: set[str], path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(6, 4.5))
    markers = dict(zip(metadata["cell"].cat.categories, "osD^v<>p"))
    colors = {"untrt": "#4C78A8", "trt": "#E45756"}
    for s, row in scores.iterrows():
        m = metadata.loc[s]
        ax.scatter(row.PC1, row.PC2, c=colors[m.dex], marker=markers[m.cell], s=70,
                   edgecolors="black" if s in flagged else "none", linewidths=2)
    for d, c in colors.items():
        ax.scatter([], [], c=c, label=f"dex={d}")
    for cell, mk in markers.items():
        ax.scatter([], [], c="gray", marker=mk, label=cell)
    ax.set_xlabel(f"PC1 ({explained[0]:.0%} of variance)")
    ax.set_ylabel(f"PC2 ({explained[1]:.0%} of variance)")
    ax.set_title("PCA of log-CPM (QC only)")
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def _plot_correlation(corr: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(corr, cmap="viridis")
    ax.set_xticks(range(len(corr)), corr.columns, rotation=90, fontsize=8)
    ax.set_yticks(range(len(corr)), corr.index, fontsize=8)
    for i in range(len(corr)):
        for j in range(len(corr)):
            ax.text(j, i, f"{corr.iat[i, j]:.3f}", ha="center", va="center", fontsize=6,
                    color="white")
    fig.colorbar(im, ax=ax, label="Pearson r (log-CPM)")
    ax.set_title("Sample-to-sample correlation")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def run_qc(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    out_dir: Path = config.RESULTS / "qc",
    rules: dict[str, Any] = config.QC_RULES,
    low_count: dict[str, int] = config.LOW_COUNT_FILTER,
) -> dict[str, Any]:
    """Run QC on raw counts and decide whether interpretation may proceed.

    Returns {passed, flags, summary, metrics, artifacts}. `passed` is False if
    any sample trips a gate rule; each flag says which rule, the value, and the
    threshold so the decision is explainable.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    lib = counts.sum()
    keep = low_count_mask(counts, **low_count)
    values = log_cpm(counts.loc[keep])
    scores, explained = pca(values, rules["pca_top_genes"])

    lib_flags = library_size_flags(lib, rules["min_library_fraction"])
    pca_flags, dist = pca_outlier_flags(scores, metadata["dex"], rules["pca_outlier_mad"])
    flags = lib_flags + pca_flags
    sep = separation_report(scores, explained, metadata, rules["separation_r2"])
    corr = values.corr()

    artifacts = {
        "library_sizes": out_dir / "library_sizes.png",
        "pca": out_dir / "pca.png",
        "correlation": out_dir / "correlation.png",
    }
    _plot_library_sizes(lib, rules["min_library_fraction"] * lib.median(), artifacts["library_sizes"])
    _plot_pca(scores, explained, metadata, {f["sample"] for f in flags}, artifacts["pca"])
    labels = metadata["title"] if "title" in metadata else pd.Series(corr.index, index=corr.index)
    _plot_correlation(corr.rename(index=labels, columns=labels), artifacts["correlation"])

    passed = not flags
    summary = (
        f"QC {'PASSED' if passed else 'FAILED'}: {len(counts.columns)} samples, "
        f"{int(keep.sum()):,} of {len(counts):,} genes pass the low-count filter "
        f"(>= {low_count['min_count']} reads in >= {low_count['min_samples']} samples). "
        f"Library sizes {lib.min() / 1e6:.1f}M-{lib.max() / 1e6:.1f}M. "
        f"PC1 explains {explained[0]:.0%} of variance (R^2 dex={sep['PC1']['r2_dex']:.2f}, "
        f"cell={sep['PC1']['r2_cell']:.2f}); PC2 {explained[1]:.0%} "
        f"(R^2 dex={sep['PC2']['r2_dex']:.2f}, cell={sep['PC2']['r2_cell']:.2f}). "
        f"Samples separate by: {sep['separates_by']}. "
        f"Sample closest to the PCA outlier cutoff: {dist.idxmax()} "
        f"({dist.max():.1f} vs cutoff {dist.attrs['cutoff']:.1f})."
    )
    if flags:
        summary += " Flags: " + " ".join(f["message"] for f in flags)

    return {
        "passed": passed,
        "flags": flags,
        "summary": summary,
        "metrics": {
            "library_sizes": {s: int(n) for s, n in lib.items()},
            "genes_total": int(len(counts)),
            "genes_passing_filter": int(keep.sum()),
            "pca_scores": scores[["PC1", "PC2"]].round(4).to_dict(orient="index"),
            "pca_centroid_distance": dist.round(4).to_dict(),
            "pca_outlier_cutoff": round(dist.attrs["cutoff"], 4),
            "separation": sep,
            "min_sample_correlation": float(corr.where(~np.eye(len(corr), dtype=bool)).min().min()),
        },
        "artifacts": {k: str(v) for k, v in artifacts.items()},
    }
