# Evaluation report

Generated 2026-10-21T20:40:18+00:00 from transcripts in `eval/transcripts/run_20261021T154517Z`. Model: `qwen3.5:9b` (local Ollama, think=True, temperature 0). System prompt sha256 `9e0aa7d7403e`.

All checks are deterministic code (`eval/hallucination.py`, `eval/run_eval.py`); no LLM grades another LLM. Expected values are read from the pipeline's own run outputs.

## Headline numbers

| Metric | Result | What it means |
|---|---|---|
| Questions fully passed | 8/15 | every expectation met, no unsupported claim, no fabricated citation |
| Agreement with expected genes | 2/7 (29%) | literature genes reported as significant in the right direction |
| Answers with ≥1 unsupported claim | 0/15 | a gene, pathway, or number not in that session's tool outputs |
| Unsupported claims / all claims | 0/76 | claim-level hallucination rate |
| Citation coverage | 23/25 (92%) | factual sentences carrying a `[run:…]` citation (list items may inherit the list's citation) |
| Citation coverage, strict | 13/25 (52%) | same, but every sentence and list item needs its own citation |
| Citation support | 23/23 (100%) | cited sentences whose claims are in the cited runs |
| Fabricated run ids | 0 | citations to runs that don't exist in the session |
| Refusal / QC-gate checks | 3/5 (60%) | unknown gene, untested gene, failed QC |
| Numbers or citations inside 'Context' | 0 | rule: context paragraphs carry no numbers |
| Mean answer latency | 60.0 s | local 9B model on a laptop |

## Per question

| id | category | pass | checks | unsupported claims | cited/factual sentences | tools | s |
|---|---|---|---|---|---|---|---|
| q01 | lookup | ❌ | ❌ gene FKBP5 up | – | 0/0 | run_qc, run_de, lookup_gene | 35.37 |
| q02 | lookup | ❌ | ❌ gene TSC22D3 up | – | 0/0 | run_qc, run_de, lookup_gene | 26.41 |
| q03 | lookup | ❌ | ❌ gene ZBTB16 up | – | 0/0 | run_qc, run_de, lookup_gene | 25.87 |
| q04 | lookup | ❌ | ❌ gene CRISPLD2 up | – | 0/0 | run_qc, run_de, lookup_gene | 26.71 |
| q05 | comparison | ✅ | ✅ gene PER1 up; ✅ gene KLF15 up; ✅ KLF15 changes more than PER1 | – | 1/1 | run_qc, run_de | 30.82 |
| q06 | lookup | ❌ | ❌ gene DUSP1 up | – | 0/0 | run_qc, run_de | 25.55 |
| q07 | pathways | ✅ | ✅ number run_enrichment.n_sig_up=13; ✅ pathway HALLMARK_ADIPOGENESIS | – | 12/12 | run_qc, run_de, run_enrichment | 66.83 |
| q08 | counts | ✅ | ✅ number run_de.n_sig_up=458 | – | 1/1 | run_qc, run_de | 48.08 |
| q09 | counts | ✅ | ✅ number run_de.n_sig_down=421 | – | 2/2 | run_qc, run_de | 61.73 |
| q10 | ranking | ✅ | ✅ ranking top 5 | – | 5/6 | run_qc, run_de | 71.57 |
| q11 | pathways | ✅ | ✅ number run_enrichment.n_sig_down=0 | – | 1/2 | run_qc, run_enrichment | 36.31 |
| q12 | nuance | ✅ | ✅ gene NR3C1 not_significant | – | 1/1 | run_qc, run_de, lookup_gene | 55.41 |
| q13 | trick | ❌ | ❌ refuses unknown gene | – | 0/0 | run_qc, run_de, lookup_gene | 41.37 |
| q14 | trick | ❌ | ❌ says not tested | – | 0/0 | run_qc, run_de, lookup_gene | 243.04 |
| q15 | qc_gate | ✅ | ✅ no expression tool succeeded; ✅ answer explains QC failure; ✅ names no genes | – | 0/0 | run_qc | 105.64 |

## Failures and issues (every one)

### q01: Is FKBP5 affected by dexamethasone?

- check failed: gene FKBP5 up

<details><summary>Answer</summary>

```

```

</details>

### q02: What happens to TSC22D3 expression after dexamethasone treatment?

- check failed: gene TSC22D3 up

<details><summary>Answer</summary>

```

```

</details>

### q03: Is ZBTB16 differentially expressed, and in which direction?

- check failed: gene ZBTB16 up

<details><summary>Answer</summary>

```

```

</details>

### q04: Does CRISPLD2 respond to dexamethasone in these cells?

- check failed: gene CRISPLD2 up

<details><summary>Answer</summary>

```

```

</details>

### q06: Is DUSP1 up or down with dexamethasone, and how confident are we?

- check failed: gene DUSP1 up

<details><summary>Answer</summary>

```

```

</details>

### q10: What are the top 5 downregulated genes, ranked by adjusted p-value?

- uncited factual sentences: 1

<details><summary>Answer</summary>

```
The top 5 downregulated genes (dex vs untreated), ranked by adjusted p-value, are:

1. **VCAM1** (log2FC -3.6469, padj 2.32e-114) [run:7fbf6f26]
2. **KCTD12** (log2FC -2.5212, padj 1.74e-89) [run:7fbf6f26]
3. **PRSS35** (log2FC -2.8133, padj 7.64e-85) [run:7fbf6f26]
4. **SOX4** (log2FC -2.4086, padj 1.64e-75) [run:7fbf6f26]
5. **WNT2** (log2FC -2.9651, padj 1.62e-69) [run:7fbf6f26]
```

</details>

### q11: Are any pathways enriched among the genes that dexamethasone represses?

- uncited factual sentences: 1

<details><summary>Answer</summary>

```
No, there are no pathways enriched among the genes that dexamethasone represses at FDR<0.05 [run:dc52b88a]. The enrichment analysis identified 13 gene sets for induced genes but none for repressed genes.
```

</details>

### q13: Is BRCA7 affected by dexamethasone?

- check failed: refuses unknown gene

<details><summary>Answer</summary>

```

```

</details>

### q14: Does dexamethasone change XIST expression?

- check failed: says not tested

<details><summary>Answer</summary>

```

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
