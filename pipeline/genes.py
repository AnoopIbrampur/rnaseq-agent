"""Gene lookup tool: one gene's DE result, or exactly why there isn't one.

A lookup can end in one of five statuses, so the agent never has to guess:

* tested: the gene has a DE row (significant or not).
* low_count_filtered: present in the count matrix but too few reads to test.
* not_in_count_matrix: a real NCBI gene that NCBI's counts don't include.
* alias: the query is an alternative name; we resolve it and report under
  the official symbol.
* unknown_symbol: not a human gene symbol or alias we know about.
"""

from __future__ import annotations

from typing import Any

import pandas as pd

from pipeline.enrichment import read_gmt


def _resolve(query: str, gene_map: pd.DataFrame) -> tuple[str | None, str | None]:
    """Return (official symbol, alias it was resolved from) for a query."""
    q = query.strip().upper()
    current = gene_map[gene_map["symbol_source"] == "gene_info"]
    if (current["symbol"].str.upper() == q).any():
        return current.loc[current["symbol"].str.upper() == q, "symbol"].iloc[0], None
    if (gene_map["symbol"].str.upper() == q).any():
        return gene_map.loc[gene_map["symbol"].str.upper() == q, "symbol"].iloc[0], None
    syn = current["synonyms"].fillna("").str.upper().str.split("|")
    hits = current.loc[syn.map(lambda s: q in s), "symbol"].unique().tolist()
    if len(hits) == 1:
        return hits[0], query
    return None, ", ".join(hits) if hits else None


def lookup_gene(
    query: str,
    de_table: pd.DataFrame,
    filtered: pd.DataFrame,
    gsea_table: pd.DataFrame,
    gene_sets: dict[str, list[str]],
    gene_map: pd.DataFrame,
    fdr_threshold: float,
    alpha: float = 0.05,
    lfc_threshold: float = 1.0,
    min_samples: int = 4,
) -> dict[str, Any]:
    """Look up a gene by symbol (or alias) in this session's DE and GSEA results."""
    symbol, alias_info = _resolve(query, gene_map)
    if symbol is None:
        msg = f"'{query}' is not a known human gene symbol."
        if alias_info:
            msg += f" It is an alias of several genes ({alias_info}); ask about one of them."
        return {"query": query, "status": "unknown_symbol", "summary": msg}

    base: dict[str, Any] = {"query": query, "symbol": symbol}
    if alias_info:
        base["resolved_from_alias"] = alias_info

    rows = de_table[de_table["symbol"] == symbol]
    if rows.empty:
        frows = filtered[filtered["symbol"] == symbol]
        if not frows.empty:
            r = frows.sort_values("max_count").iloc[-1]
            return {**base, "status": "low_count_filtered",
                    "max_count": int(r["max_count"]),
                    "n_samples_passing": int(r["n_samples_passing"]),
                    "summary": (f"{symbol} was not tested: it was removed by the low-count "
                                f"filter. Only {int(r['n_samples_passing'])} of 8 samples had "
                                f"enough reads (highest count in any sample: "
                                f"{int(r['max_count'])}); the filter requires "
                                f"{min_samples} samples.")}
        return {**base, "status": "not_in_count_matrix",
                "summary": f"{symbol} is not in this dataset's count matrix, so it was not tested."}

    current = rows[rows["symbol_source"] == "gene_info"]
    row = (current if not current.empty else rows).sort_values("baseMean").iloc[-1]
    member_of = sorted(s for s, genes in gene_sets.items() if symbol in genes)
    g = gsea_table.set_index("term")
    sets = []
    for s in member_of:
        if s not in g.index:
            continue
        r = g.loc[s]
        sets.append({
            "term": s, "nes": round(float(r["nes"]), 4), "fdr": float(r["fdr"]),
            "significant": bool(r["fdr"] < fdr_threshold),
            "in_leading_edge": symbol in str(r["leading_edge"]).split(";"),
        })
    padj = None if pd.isna(row["padj"]) else float(row["padj"])
    result = {
        **base,
        "status": "tested",
        "gene_id": int(row["gene_id"]),
        "baseMean": round(float(row["baseMean"]), 2),
        "log2FC": round(float(row["log2FC"]), 4),
        "lfcSE": round(float(row["lfcSE"]), 4),
        "stat": round(float(row["stat"]), 4),
        "padj": padj,
        "significant": bool(row["significant"]),
        "direction": row["direction"],
        "gene_sets": sets,
        "other_rows_with_same_symbol": int(len(rows) - 1),
    }
    word = "higher" if row["log2FC"] > 0 else "lower"
    if row["significant"]:
        verdict = f"significantly {word} with dex (padj<{alpha}, |log2FC|>{lfc_threshold:g})"
    elif padj is not None and padj < alpha:
        verdict = (f"statistically {word} with dex (padj<{alpha}) but the change is below the "
                   f"|log2FC|>{lfc_threshold:g} threshold, so not counted as significant")
    else:
        verdict = f"no significant change (padj>={alpha})"
    result["significance_detail"] = verdict
    sig_sets = [s["term"] for s in sets if s["significant"]]
    result["summary"] = (
        f"{symbol}: log2FC={result['log2FC']:+.2f} (shrunken), padj={padj:.3g}, {verdict}. "
        + (f"Member of significant gene sets: {', '.join(sig_sets)}." if sig_sets
           else "Not a member of any significant gene set.")
    ) if padj is not None else f"{symbol}: tested but no adjusted p-value ({row['note']})."
    return result


__all__ = ["lookup_gene", "read_gmt"]
