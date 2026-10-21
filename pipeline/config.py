"""Central configuration: paths, expected checksums, analysis parameters.

Everything that can change a result lives here so it can be recorded in
provenance records and hashed into cache keys.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA_RAW = ROOT / "data" / "raw"
DATA_PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
RUNS = ROOT / "runs"

# Raw inputs, downloaded once and checksummed (see data/SOURCES.md).
COUNTS_FILE = DATA_RAW / "GSE52778_raw_counts_GRCh38.p13_NCBI.tsv.gz"
SERIES_MATRIX_FILE = DATA_RAW / "GSE52778_series_matrix.txt.gz"
GENE_INFO_FILE = DATA_RAW / "Homo_sapiens.gene_info.gz"
AUTHORS_DE_FILE = DATA_RAW / "GSE52778_Dex_vs_Untreated_gene_exp.diff.gz"
# Human rows of NCBI gene_history: maps retired GeneIDs to current ones.
GENE_HISTORY_FILE = DATA_RAW / "gene_history_human.tsv.gz"

EXPECTED_SHA256: dict[str, str] = {
    COUNTS_FILE.name: "8f5d0adafcb5d943338e196b99aa87dd8ca2515927450627b7443328a3956603",
    SERIES_MATRIX_FILE.name: "67db46e962f795b9dd61fa6b4279f2678c6b2b9bf1df035fbc22fcf81b5e879e",
    GENE_INFO_FILE.name: "0bab97eb8e0373c91d092b95259b8389ffd6a57a5c9f146bfe540968ab39542c",
    AUTHORS_DE_FILE.name: "6c6ca5ce509e98ea9d5fb9aef698b2a0fbf4c778997c7ceaf95a42cbf15eb63f",
    GENE_HISTORY_FILE.name: "a7a815e4b3517f5f0f7193cd461d2d2ff5d505ed8295a910e2325e31284fa542",
}

# Derived files (rebuilt from raw inputs, safe to delete).
GENE_MAP_FILE = DATA_PROCESSED / "gene_map.csv"

# GEO holds 16 samples: 4 donor cell lines x {Untreated, Dexamethasone,
# Albuterol, Albuterol_Dexamethasone}. This project studies one contrast,
# dex vs untreated, so we keep only those 8 samples (the same 8 as the
# Bioconductor `airway` package).
KEEP_TREATMENTS = {"Untreated": "untrt", "Dexamethasone": "trt"}

# Some retired GeneIDs in the count matrix were merged by NCBI into a gene
# that also has its own row (e.g. LOC112268437 -> AAK1), so ~320 symbols
# appear on 2+ rows. DE runs per GeneID on NCBI's counts unchanged. For
# symbol-keyed steps (GSEA ranking, gene lookup) we keep one row per symbol:
# the current gene_info GeneID if present, else the highest-baseMean row.
SYMBOL_COLLAPSE_RULE = "prefer_current_geneid_then_max_basemean"

SEED = 42
