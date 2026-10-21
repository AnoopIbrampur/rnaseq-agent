# Evaluation report

Generated 2026-10-21T20:51:40+00:00 from transcripts in `eval/transcripts/run_20261021T204031Z`. Model: `qwen3.5:9b` (local Ollama, think=True, temperature 0). System prompt sha256 `9e0aa7d7403e`.

All checks are deterministic code (`eval/hallucination.py`, `eval/run_eval.py`); no LLM grades another LLM. Expected values are read from the pipeline's own run outputs.

## Headline numbers

| Metric | Result | What it means |
|---|---|---|
| Questions fully passed | 14/15 | every expectation met, no unsupported claim, no fabricated citation |
| Agreement with expected genes | 7/7 (100%) | literature genes reported as significant in the right direction |
| Answers with ≥1 unsupported claim | 1/15 | a gene, pathway, or number not in that session's tool outputs |
| Unsupported claims / all claims | 2/120 | claim-level hallucination rate |
| Citation coverage | 29/34 (85%) | factual sentences carrying a `[run:…]` citation (list items may inherit the list's citation) |
| Citation coverage, strict | 19/34 (56%) | same, but every sentence and list item needs its own citation |
| Citation support | 29/29 (100%) | cited sentences whose claims are in the cited runs |
| Fabricated run ids | 0 | citations to runs that don't exist in the session |
| Refusal / QC-gate checks | 5/5 (100%) | unknown gene, untested gene, failed QC |
| Numbers or citations inside 'Context' | 0 | rule: context paragraphs carry no numbers |
| Mean answer latency | 44.5 s | local 9B model on a laptop |

## Per question

| id | category | pass | checks | unsupported claims | cited/factual sentences | tools | s |
|---|---|---|---|---|---|---|---|
| q01 | lookup | ✅ | ✅ gene FKBP5 up | – | 2/2 | run_qc, run_de, lookup_gene | 58.36 |
| q02 | lookup | ✅ | ✅ gene TSC22D3 up | – | 1/1 | run_qc, run_de, lookup_gene | 48.15 |
| q03 | lookup | ✅ | ✅ gene ZBTB16 up | – | 1/1 | run_qc, run_de, lookup_gene | 44.77 |
| q04 | lookup | ✅ | ✅ gene CRISPLD2 up | – | 1/1 | run_qc, run_de, lookup_gene | 50.96 |
| q05 | comparison | ✅ | ✅ gene PER1 up; ✅ gene KLF15 up; ✅ KLF15 changes more than PER1 | – | 1/1 | run_qc, run_de | 27.29 |
| q06 | lookup | ✅ | ✅ gene DUSP1 up | – | 1/1 | run_qc, run_de | 46.03 |
| q07 | pathways | ✅ | ✅ number run_enrichment.n_sig_up=13; ✅ pathway HALLMARK_ADIPOGENESIS | – | 12/12 | run_qc, run_de, run_enrichment | 53.42 |
| q08 | counts | ✅ | ✅ number run_de.n_sig_up=458 | – | 1/1 | run_qc, run_de | 24.54 |
| q09 | counts | ✅ | ✅ number run_de.n_sig_down=421 | – | 1/2 | run_qc, run_de | 57.25 |
| q10 | ranking | ✅ | ✅ ranking top 5 | – | 5/6 | run_qc, run_de | 66.98 |
| q11 | pathways | ✅ | ✅ number run_enrichment.n_sig_down=0 | – | 1/2 | run_qc, run_enrichment | 28.67 |
| q12 | nuance | ✅ | ✅ gene NR3C1 not_significant | – | 1/1 | run_qc, run_de, lookup_gene | 39.44 |
| q13 | trick | ❌ | ✅ refuses unknown gene | BRCA1, BRCA2 | 0/1 | run_qc, run_de, lookup_gene | 50.79 |
| q14 | trick | ✅ | ✅ says not tested | – | 1/2 | run_qc, run_de, lookup_gene | 49.61 |
| q15 | qc_gate | ✅ | ✅ no expression tool succeeded; ✅ answer explains QC failure; ✅ names no genes | – | 0/0 | run_qc | 20.92 |

## Failures and issues (every one)

### q09: How many genes go down with dexamethasone?

- uncited factual sentences: 1

<details><summary>Answer</summary>

```
421 genes go down with dexamethasone [run:981c3a7d]. The top 10 include VCAM1 (log2FC -3.65, padj 2.32e-114), KCTD12 (log2FC -2.52, padj 1.74e-89), PRSS35 (log2FC -2.81, padj 7.64e-85), SOX4 (log2FC -2.41, padj 1.64e-75), WNT2 (log2FC -2.97, padj 1.62e-69), CCKAR (log2FC -2.16, padj 2.17e-69), ADAM12 (log2FC -1.82, padj 1.11e-68), ADAMTS14 (log2FC -2.47, padj 2.61e-62), CXCL12 (log2FC -1.86, padj 3.58e-62), and FER1L6 (log2FC -3.07, padj 5.1e-60).
```

</details>

### q10: What are the top 5 downregulated genes, ranked by adjusted p-value?

- uncited factual sentences: 1

<details><summary>Answer</summary>

```
The top 5 downregulated genes (ranked by adjusted p-value) are:

1. **VCAM1** (log2FC -3.65, padj 2.32e-114) [run:0901ef80]
2. **KCTD12** (log2FC -2.52, padj 1.74e-89) [run:0901ef80]
3. **PRSS35** (log2FC -2.81, padj 7.64e-85) [run:0901ef80]
4. **SOX4** (log2FC -2.41, padj 1.64e-75) [run:0901ef80]
5. **WNT2** (log2FC -2.97, padj 1.62e-69) [run:0901ef80]
```

</details>

### q11: Are any pathways enriched among the genes that dexamethasone represses?

- uncited factual sentences: 1

<details><summary>Answer</summary>

```
No, there are no pathways enriched among the genes that dexamethasone represses at FDR<0.05 [run:b2938dd5]. The analysis identified 13 enriched gene sets for dexamethasone-induced genes, but none for repressed genes.
```

</details>

### q13: Is BRCA7 affected by dexamethasone?

- unsupported genes: BRCA1, BRCA2
- uncited factual sentences: 1

<details><summary>Answer</summary>

```
BRCA7 is not a known human gene symbol in this analysis [run:df223846]. It was not tested for differential expression. If you meant a different gene (e.g., BRCA1 or BRCA2), please let me know and I can look those up.
```

</details>

### q14: Does dexamethasone change XIST expression?

- uncited factual sentences: 1

<details><summary>Answer</summary>

```
XIST was not tested in this analysis: it was removed by the low-count filter because only 2 of 8 samples had enough reads (the filter requires ≥4 samples) [run:b16cb5e6]. Therefore, no differential expression result is available for XIST with dexamethasone.
```

</details>


## Pipeline sanity checks (no LLM involved)

### Expected genes vs this pipeline's DE results

Sanity check, not ground truth: see `eval/expected.yaml` for where the list comes from.

| gene | expected | observed | log2FC (shrunken) | padj | ok |
|---|---|---|---|---|---|
| FKBP5 | up | up | 3.962 | 2.1e-26 | ✅ |
| TSC22D3 | up | up | 3.134 | 2.1e-19 | ✅ |
| PER1 | up | up | 3.207 | 6.3e-99 | ✅ |
| KLF15 | up | up | 4.395 | 1.2e-113 | ✅ |
| DUSP1 | up | up | 2.946 | 1.1e-139 | ✅ |
| ZBTB16 | up | up | 7.376 | 1.1e-58 | ✅ |
| CRISPLD2 | up | up | 2.61 | 3.2e-61 | ✅ |

### Agreement with the study authors' own analysis

The authors' Cuffdiff dex-vs-untreated table from GEO (`GSE52778_Dex_vs_Untreated_gene_exp.diff.gz`) is an independent pipeline (different aligner, annotation, statistical model, and no paired design), so perfect overlap is not expected.

- Genes matched by symbol: 11,379
- Significant here: 876 (padj<0.05, |log2FC|>1); significant in Cuffdiff: 316 (its own q<0.05 call)
- Significant in both: 280; same direction: 280
- Jaccard overlap of significant sets: 0.307
- Spearman correlation of log2 fold changes (unshrunken here): 0.959

## Method notes and known blind spots

- **Hallucination** = a gene symbol, MSigDB set name, or number in the findings part of an answer that does not appear in anything the agent was shown in that session (tool results, the question, the system prompt). Numbers match with rounding tolerance (half a unit in the last written decimal; 5% relative for scientific notation).
- **Not detected:** pathway names written in plain English, spelled-out numbers ("thirteen"), and wrong *relationships* between correct facts (e.g. right number, wrong gene). The ranking check covers one such case explicitly.
- **Direction claims** are read with keyword rules (higher/up/induced vs lower/down/repressed, plus 'not significant' phrases). Ambiguous sentences are scored `unclear`, which counts as wrong.
- **Citation support** re-checks each cited sentence against only the runs it cites (their results plus the relevant rows of their output tables).
- One run of 15 questions at temperature 0 is a small sample; rerunning can change answers. See `eval/CHANGELOG.md` for any prompt changes made after looking at results.
