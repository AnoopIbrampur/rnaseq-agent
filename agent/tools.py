"""Session-level tool implementations: what the agent (and CLI) actually call.

A Session fixes the data source and remembers this session's runs, so
downstream tools reuse the upstream run (e.g. enrichment cites the DE run it
ranked) instead of silently recomputing. Every call goes through
provenance.execute(), so every result has a run record.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from pipeline import config, provenance


@dataclass
class DataSource:
    """Which count matrix and metadata a session analyzes."""

    counts: Path = config.COUNTS_FILE
    metadata: Path = config.SERIES_MATRIX_FILE


@dataclass
class Session:
    """One analysis session: a data source plus the run records produced so far."""

    source: DataSource = field(default_factory=DataSource)
    runs: dict[str, dict[str, Any]] = field(default_factory=dict)
    latest: dict[str, str] = field(default_factory=dict)

    def _record(self, rec: dict[str, Any]) -> dict[str, Any]:
        self.runs[rec["run_id"]] = rec
        self.latest[rec["tool"]] = rec["run_id"]
        return rec

    def _latest(self, tool: str) -> dict[str, Any] | None:
        rid = self.latest.get(tool)
        return self.runs[rid] if rid else None

    def run_qc(self) -> dict[str, Any]:
        """QC with gate on this session's data."""
        return self._record(provenance.execute(
            "run_qc",
            {"counts": self.source.counts, "metadata": self.source.metadata},
            {"qc_rules": config.QC_RULES, "low_count": config.LOW_COUNT_FILTER},
        ))

    def run_de(self) -> dict[str, Any]:
        """Differential expression, dex vs untreated, design ~cell + dex."""
        if (rec := self._latest("run_de")) is not None:
            return rec
        return self._record(provenance.execute(
            "run_de",
            {"counts": self.source.counts, "metadata": self.source.metadata,
             "gene_info": config.GENE_INFO_FILE, "gene_history": config.GENE_HISTORY_FILE},
            {"de": config.DE_PARAMS, "low_count": config.LOW_COUNT_FILTER},
        ))

    def run_enrichment(self) -> dict[str, Any]:
        """Prerank GSEA on this session's DE results."""
        if (rec := self._latest("run_enrichment")) is not None:
            return rec
        de = self.run_de()
        return self._record(provenance.execute(
            "run_enrichment",
            {"de_table": provenance.absolute(de["outputs"]["table"]),
             "gene_sets": config.GENE_SETS_FILE},
            config.GSEA_PARAMS,
            upstream_runs=[de["run_id"]],
        ))

    def lookup_gene(self, symbol: str) -> dict[str, Any]:
        """One gene's DE result and gene-set membership, or why it has none."""
        de, gsea = self.run_de(), self.run_enrichment()
        return self._record(provenance.execute(
            "lookup_gene",
            {"de_table": provenance.absolute(de["outputs"]["table"]),
             "filtered_genes": provenance.absolute(de["outputs"]["filtered_genes"]),
             "gsea_table": provenance.absolute(gsea["outputs"]["table"]),
             "gene_sets": config.GENE_SETS_FILE,
             "gene_info": config.GENE_INFO_FILE, "gene_history": config.GENE_HISTORY_FILE},
            {"query": symbol, "fdr_threshold": config.GSEA_PARAMS["fdr_threshold"],
             "alpha": config.DE_PARAMS["alpha"],
             "lfc_threshold": config.DE_PARAMS["lfc_threshold"],
             "min_samples": config.LOW_COUNT_FILTER["min_samples"]},
            upstream_runs=[de["run_id"], gsea["run_id"]],
        ))
