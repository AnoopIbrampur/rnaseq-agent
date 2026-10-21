"""Build notebooks/01_walkthrough.ipynb from cells defined here, then execute it.

    .venv/bin/python notebooks/build_walkthrough.py

Keeping the notebook as code makes it reviewable in diffs; the executed .ipynb
(with outputs) is what a reader opens.
"""

from pathlib import Path

import nbformat as nbf
from nbconvert.preprocessors import ExecutePreprocessor

HERE = Path(__file__).resolve().parent
md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell

cells = [
    md("""# Airway RNA-seq: the analysis without the agent

This notebook runs the same pipeline functions the agent calls (`pipeline/`), step by step, so a computational biologist can check the science directly. No language model is involved.

**Experiment.** Human airway smooth muscle cells from 4 donors, each untreated or treated with dexamethasone (a glucocorticoid). GEO [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778), Himes et al. 2014. We use NCBI's reprocessed raw counts (GRCh38.p13) and keep the 8 untreated/dex samples (GEO also has albuterol arms).

**Sections:** 1 data · 2 QC · 3 differential expression · 4 why the paired design matters · 5 pathway enrichment · 6 checks against the literature genes and against the authors' own analysis · 7 edge cases the agent must handle"""),
    code("""import sys, warnings, logging
from pathlib import Path
ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
logging.disable(logging.WARNING)  # gseapy logs tie warnings; ties are broken deterministically

import pandas as pd, numpy as np
from IPython.display import Image, display
from pipeline import config, data, qc, de, enrichment
pd.set_option("display.width", 140, "display.max_columns", 12)
OUT = config.RESULTS / "notebook"
print("pydeseq2", __import__("pydeseq2").__version__, "| gseapy", __import__("gseapy").__version__)"""),

    md("""## 1. Data

Raw files are checksummed on every load against `data/SOURCES.md`. The count matrix holds **raw integer read counts**: DESeq2 models counts directly, so we never feed it normalized values (FPKM/TPM)."""),
    code("""ds = data.load_dataset()
print(ds.counts.shape, "genes x samples")
ds.metadata"""),
    code("""lib = ds.counts.sum()
pd.DataFrame({"reads (M)": (lib / 1e6).round(1), "cell": ds.metadata.cell, "dex": ds.metadata.dex})"""),

    md("""## 2. Quality control

Three explainable gate rules (in `config.QC_RULES`): library size ≥ 50% of the median; no PCA outlier beyond 3 robust SDs from its treatment group's centroid; every sample's best-matching other sample correlates at r ≥ 0.90. PCA and correlations use log-CPM of genes passing the low-count filter, for visualization only."""),
    code("""q = qc.run_qc(ds.counts, ds.metadata, out_dir=OUT / "qc")
print(q["summary"])
q["flags"]"""),
    code("""display(Image(q["artifacts"]["pca"], width=520))"""),
    md("""PC1 is treatment (dex explains ~95% of its variance), PC2 is donor (cell line explains ~99%). Treatment is the dominant signal, and the donor effect is large: that is what the `cell` term in the design is for.

**The gate on a broken sample.** Scrambling one sample's counts across genes (simulating a mislabeled or degraded library) trips both the PCA and correlation rules:"""),
    code("""bad = data.corrupt_counts(ds.counts, "scrambled_sample", "GSM1275867")
qb = qc.run_qc(bad, ds.metadata, out_dir=OUT / "qc_corrupted")
print(qb["passed"]); [f["message"] for f in qb["flags"]]"""),

    md("""## 3. Differential expression

PyDESeq2 with design `~ cell + dex`, Wald test for dex (treated vs untreated), Benjamini–Hochberg adjusted p-values, and LFC shrinkage. "Significant" means padj < 0.05 and |shrunken log2FC| > 1. The reference level is set explicitly to `untrt`: with string labels PyDESeq2 would pick `trt` alphabetically and the shrunken fold changes would flip sign."""),
    code("""d = de.run_de(ds.counts, ds.metadata, ds.gene_map, out_dir=OUT / "de")
print(d["summary"]); print(d["warnings"])
table = pd.read_csv(d["artifacts"]["table"])
table[table.significant].head(20)[["symbol", "baseMean", "log2FC", "log2FC_mle", "stat", "padj", "direction"]]"""),
    code("""display(Image(d["artifacts"]["volcano"], width=560))"""),
    md("""The small-sample warning is expected: 8 samples minus 5 model coefficients leaves 3 residual degrees of freedom. Shrinkage barely moves the estimates here (compare `log2FC` with `log2FC_mle`), a sign the counts are deep and consistent."""),

    md("""## 4. Why the paired design matters

Same data, design `~ dex` (ignoring which donor each sample came from). Donor differences then count as noise, inflating the variance of every gene."""),
    code("""unpaired = {**config.DE_PARAMS, "design": "~dex"}
d0 = de.run_de(ds.counts, ds.metadata, ds.gene_map, out_dir=OUT / "de_unpaired", params=unpaired)
t0 = pd.read_csv(d0["artifacts"]["table"]).drop_duplicates("gene_id").set_index("gene_id")
t1 = table.drop_duplicates("gene_id").set_index("gene_id")
pd.DataFrame({
    "~cell + dex (paired)": [d["n_sig_up"], d["n_sig_down"], round(t1.stat.abs().median(), 2)],
    "~dex (unpaired)": [d0["n_sig_up"], d0["n_sig_down"], round(t0.stat.abs().median(), 2)],
}, index=["significant up", "significant down", "median |Wald stat|"])"""),

    md("""## 5. Pathway enrichment

Prerank GSEA (gseapy) of all tested genes ranked by Wald statistic against MSigDB Hallmark v2026.2 (50 sets, local file). Unlike over-representation analysis, it needs no significance cutoff and can detect many genes shifting modestly together. Duplicate symbols (from retired NCBI GeneIDs) are collapsed to one row per symbol, preferring the current GeneID."""),
    code("""g = enrichment.run_enrichment(table, out_dir=OUT / "gsea")
print(g["summary"])
gt = pd.read_csv(g["artifacts"]["table"])
gt.head(15)[["term", "nes", "fdr", "set_size", "leading_edge"]].assign(
    leading_edge=lambda x: x.leading_edge.str.split(";").str[:6].str.join(", "))"""),
    code("""display(Image(g["artifacts"]["plot"], width=560))"""),
    md("""**Read the leading edge, not just the name.** `TNFA_SIGNALING_VIA_NFKB` is enriched among *dex-induced* genes, which sounds backwards for an anti-inflammatory drug. Its leading edge is DUSP1, PER1, KLF9 and similar genes. DUSP1 and PER1 are named as well-known glucocorticoid-responsive genes in this dataset's paper (Himes et al. 2014), and KLF9 is glucocorticoid-induced with GR binding sites upstream (Mostafa et al. 2021, *J Biol Chem*). *(Interpretation, not computed here: glucocorticoid targets that sit in the set drive its score.)* This is why the agent cites leading-edge genes."""),

    md("""## 6. Checks

### 6a. Literature genes (sanity check, not ground truth)
Glucocorticoid-responsive genes with literature support: FKBP5, TSC22D3, PER1, KLF15, DUSP1, CRISPLD2 from this dataset's paper (Himes et al. 2014, PMID 24926665); ZBTB16 from an independent airway smooth muscle study (Kan et al. 2019, PMID 30694689). Sources per gene are in `eval/expected.yaml`."""),
    code("""import yaml
exp = yaml.safe_load(open(ROOT / "eval/expected.yaml"))["genes"]
t = table[table.symbol_source == "gene_info"].set_index("symbol")
t.loc[list(exp), ["log2FC", "padj", "direction"]].assign(expected=list(exp.values()))"""),
    md("""### 6b. Agreement with the study authors' own results
The authors deposited a Cuffdiff dex-vs-untreated table. It is an independent pipeline (different aligner, annotation, model, and unpaired), so the overlap shouldn't be expected to be perfect. Cuffdiff's `log2(fold_change)` is log2(sample_2/sample_1) with sample_1 = Dex, so we flip its sign."""),
    code("""a = pd.read_csv(config.AUTHORS_DE_FILE, sep="\\t")
a = a[a.status == "OK"].drop_duplicates("gene").set_index("gene")
a_lfc = -pd.to_numeric(a["log2(fold_change)"], errors="coerce").replace([np.inf, -np.inf], np.nan)
ours = table.drop_duplicates("symbol").set_index("symbol")
common = ours.index.intersection(a.index)
both = [x for x in common if ours.loc[x, "significant"] and a.loc[x, "significant"] == "yes"]
print("genes in common:", len(common))
print("significant here:", int(ours.significant.sum()), "| significant in Cuffdiff:", int((a.significant == "yes").sum()))
print("significant in both:", len(both), "| same direction:", sum((ours.loc[x, "log2FC"] > 0) == (a_lfc[x] > 0) for x in both))
pair = pd.concat([ours.loc[common, "log2FC_mle"], a_lfc[common]], axis=1).dropna()
print("Spearman r of log2FC:", round(pair.corr(method="spearman").iloc[0, 1], 3))"""),

    md("""## 7. Edge cases the agent must handle

* **NR3C1** (the glucocorticoid receptor): highly significant padj but |log2FC| just under 1, so *not* counted as significant. The lookup says so explicitly.
* **XIST**: thousands of reads in one donor, near zero in the others, so it fails the low-count filter (needs ≥10 reads in ≥4 samples) and is reported as *not tested*. A female-specific gene would show this pattern. *(Context, not computed.)*"""),
    code("""print(t.loc["NR3C1", ["log2FC", "padj", "significant"]].to_dict())
ids = ds.gene_map.loc[ds.gene_map.symbol == "XIST", "gene_id"]
xist = ids[ids.isin(ds.counts.index)].iloc[0]  # XIST has >1 NCBI GeneID; take the counted one
pd.concat([ds.counts.loc[xist].rename("XIST reads"), ds.metadata.cell], axis=1)"""),

    md("""## Reproducibility

Each of these steps, when run through `provenance.execute()` (as the agent and `cli.py` do), writes a run record with input checksums, parameters, software versions and output checksums. `python cli.py reproduce <run_id>` re-runs a record from scratch and verifies the outputs match; in testing, DE and GSEA outputs reproduced byte-for-byte."""),
]

nb = nbf.v4.new_notebook(cells=cells, metadata={
    "kernelspec": {"name": "rnaseq-agent", "display_name": "rnaseq-agent (.venv)", "language": "python"}})
ExecutePreprocessor(timeout=900, kernel_name="rnaseq-agent").preprocess(nb, {"metadata": {"path": str(HERE)}})
out = HERE / "01_walkthrough.ipynb"
nbf.write(nb, out)
print(f"wrote {out}")
