# Data sources

All raw files are in `data/raw/` and committed to git (7.6 MB total), so the project runs offline.
`pipeline/data.py` checks each file's SHA-256 against `pipeline/config.py:EXPECTED_SHA256` every time it loads, and stops with an error on a mismatch.

Downloaded **2026-10-21** (UTC).

| File | What it is | Source URL | SHA-256 |
|---|---|---|---|
| `GSE52778_raw_counts_GRCh38.p13_NCBI.tsv.gz` | Raw integer read counts, 39,376 genes × 16 samples, keyed by NCBI GeneID. NCBI-generated: reads reprocessed by NCBI against GRCh38.p13, not by the study authors. | `https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE52778&format=file&file=GSE52778_raw_counts_GRCh38.p13_NCBI.tsv.gz` | `8f5d0adafcb5d943338e196b99aa87dd8ca2515927450627b7443328a3956603` |
| `GSE52778_series_matrix.txt.gz` | Sample metadata (treatment, cell line, ERCC spike-in) for all 16 samples | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52778/matrix/GSE52778_series_matrix.txt.gz` | `67db46e962f795b9dd61fa6b4279f2678c6b2b9bf1df035fbc22fcf81b5e879e` |
| `Homo_sapiens.gene_info.gz` | NCBI Gene table: GeneID → symbol, Ensembl xref, gene type. Server Last-Modified 2026-10-20 06:01 GMT (NCBI updates it daily, so this checksum pins the version). | `https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz` | `0bab97eb8e0373c91d092b95259b8389ffd6a57a5c9f146bfe540968ab39542c` |
| `gene_history_human.tsv.gz` | Human (tax_id 9606) rows of NCBI `gene_history`: maps retired GeneIDs to their replacements. **Derived:** filtered with `awk -F'\t' 'NR==1 \|\| $1=="9606"' \| gzip -n`. | `https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_history.gz` (original all-species file, Last-Modified 2026-10-20 05:53 GMT, sha256 `4ff658a71f1c1142cd67353de918439c224e6178f1e489813d34949dcc56f6ce`, not kept: 163 MB) | `a7a815e4b3517f5f0f7193cd461d2d2ff5d505ed8295a910e2325e31284fa542` |
| `GSE52778_Dex_vs_Untreated_gene_exp.diff.gz` | The study authors' own Cuffdiff dex-vs-untreated results. Not used by the pipeline; kept as an independent reference for evaluation. | `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52778/suppl/GSE52778_Dex_vs_Untreated_gene_exp.diff.gz` | `6c6ca5ce509e98ea9d5fb9aef698b2a0fbf4c778997c7ceaf95a42cbf15eb63f` |

| `h.all.v2026.2.Hs.symbols.gmt` | MSigDB **Hallmark** gene sets, release 2026.2 (human, gene symbols): 50 sets, 7,322 gene memberships. Downloaded 2026-10-21; server Last-Modified 2026-10-15 18:01 GMT. Used for prerank GSEA offline. License terms (https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp) checked by the repository owner on 2026-10-21 before publishing. Attribution: Liberzon et al., *Cell Systems* 2015 (Hallmark collection). | `https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.2.Hs/h.all.v2026.2.Hs.symbols.gmt` | `56bba1300555fb6302844fc8cdca63406eb50fa1a2472ec86468747efb3fb89b` |

Study: Himes BE et al. (2014) "RNA-Seq transcriptome profiling identifies CRISPLD2 as a glucocorticoid responsive gene that modulates cytokine function in airway smooth muscle cells." *PLoS One* 9(6):e99625. GEO GSE52778.

## Sample subset

GEO has 16 samples: 4 donor cell lines (N61311, N052611, N080611, N061011) × 4 conditions (Untreated, Dexamethasone, Albuterol, Albuterol+Dexamethasone).
We keep the **8 Untreated and Dexamethasone samples**. These are the same 8 as the Bioconductor `airway` package (`config.KEEP_TREATMENTS`).

## Gene ID → symbol mapping

`data/processed/gene_map.csv` is built by `pipeline.data.build_gene_map()`. It's derived data, not committed, and rebuilt automatically when missing. Coverage of the 39,376 count-matrix GeneIDs is 100%:

- 37,663 current GeneIDs are found directly in `gene_info`.
- 378 retired GeneIDs map to their replacement GeneID's current symbol (`symbol_source = gene_history:replaced_by_<id>`).
- 1,335 retired GeneIDs with no replacement keep their last symbol, mostly `LOC…` placeholders (`symbol_source = gene_history:discontinued`).

Because of merges, ~320 symbols appear on 2+ rows. See `config.SYMBOL_COLLAPSE_RULE` for how symbol-keyed steps choose one row.

## Known caveats

- **NCBI IDs, not Ensembl.** NCBI's pipeline uses NCBI Gene annotation. Ensembl IDs are carried as a cross-reference column where NCBI provides one.
- **ERCC spike-ins.** Some samples had ERCC synthetic RNA spiked in (`ercc_mix` column). NCBI's counts are against the human genome only, so ERCC reads aren't in the matrix. Recorded for transparency only; they have no effect on the analysis.
- **Blocked mirror.** NCBI's own `Human.GRCh38.p13.annot.tsv.gz` link returned a reCAPTCHA page to scripted downloads, so the mapping uses the FTP `gene_info` and `gene_history` files instead.
