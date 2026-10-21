# 2-minute demo

**Before the call (5 min):** warm everything up so nothing is slow live.

```bash
ollama serve                                    # if not already running
.venv/bin/python cli.py analyze                 # fills the cache, so tools answer in ~0.03 s
.venv/bin/streamlit run ui/app.py               # open http://localhost:8501
```

Ask one throwaway question in the UI so the model is loaded into memory. The first answer after startup takes the longest. Keep a terminal open in the repo, and have `eval/report.md` open in a tab.

---

## 0:00–0:20 · Ask a question

> "The idea: an agent that can't claim anything it didn't compute. This is public RNA-seq data: airway muscle cells from 4 donors, with and without dexamethasone, a steroid."

Click the example **"Is FKBP5 affected, and how confident are we?"**

> "It's a 9-billion-parameter model running locally on my laptop, no API. It can only answer by calling analysis tools, and it has to run QC first."

*(The answer takes ~30–60 s. Talk over the wait, or ask the question just before you start.)*

## 0:20–0:50 · Show the cited answer

Point at the answer:

> "Every factual sentence ends in a run id. *Higher with dexamethasone, log2 fold change +3.96, adjusted p 2e-26.* That number comes from a DESeq2 run, not from the model's memory."

Point at the caption line (`tools: run_qc(), run_de(), lookup_gene(FKBP5)`):

> "These are the tools it called. If it cites a run id that doesn't exist in this session, the UI flags it in red. If QC had failed, the interpretation tools would have refused to run, and that's enforced in code, not in the prompt."

Optional, if there's time: switch the sidebar **Data** to *Corrupted copy: one scrambled sample* and ask anything. It explains the QC flags and stops.

## 0:50–1:10 · Expand a provenance record

Click the `run:… · lookup_gene` panel:

> "This is the run record behind the claim: the input files with SHA-256 checksums, every parameter, package versions, a hash of the tool's own code, the git commit, and the output files with their checksums. It points to the upstream DE and enrichment runs, so you can follow a claim all the way back to the raw counts."

## 1:10–1:30 · Reproduce it

Copy the command at the bottom of the panel into the terminal:

```bash
.venv/bin/python cli.py reproduce <run_id>
```

> "This re-executes the run from scratch, not from cache. It checks the inputs haven't changed, then compares every output. Here it's byte-identical. The same works for the DESeq2 run itself."

*(Reproducing a lookup takes ~1 s. To reproduce the DE run instead, use its id from `cli.py runs`; that takes ~9 s.)*

## 1:30–2:00 · Show the eval table

Open `eval/report.md`, headline table:

> "I evaluate the agent with deterministic code, not another LLM. Every gene, pathway, and number in an answer must appear in that session's tool output. On 15 questions: 14 pass, 7 of 7 literature genes are reported correctly, 2 unsupported claims out of 120, and every cited sentence is actually supported by the run it cites."

> "The first run scored 7 out of 15. The model was putting its answer in its hidden reasoning and returning a blank reply. I fixed that with a retry, logged it in a changelog, and left the remaining issues documented instead of tuning the prompt until the numbers looked good."

---

## Likely follow-up questions

* **"Why a local model?"** Data never leaves the machine, it costs nothing, and the guardrails don't depend on a frontier model behaving. A bigger model would plug into the same `chat()` interface.
* **"How do you know the DE is right?"** The notebook (`notebooks/01_walkthrough.ipynb`) runs it without the agent. Fold changes correlate at Spearman 0.959 with the authors' own Cuffdiff results, and all 280 genes significant in both pipelines change in the same direction.
* **"What does the agent get wrong?"** It sometimes adds uncited restating sentences. It once suggested BRCA1/BRCA2 for a nonexistent gene. Its QC-failure message implies chat could override the gate (it can't). All of this is in `eval/CHANGELOG.md`.
* **"Why is TNFα/NF-κB *up* with an anti-inflammatory steroid?"** The leading edge is DUSP1, PER1 and KLF9, which are commonly described as glucocorticoid-induced genes that also sit in that gene set. That's why answers cite leading-edge genes and not just pathway names.

**Fallback if Ollama misbehaves live:** `bash scripts/demo.sh` shows a recorded answer from `eval/transcripts/`, a provenance record, `reproduce`, and the eval table, all without the model.
