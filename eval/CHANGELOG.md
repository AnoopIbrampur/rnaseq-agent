# Eval changelog

Any change to the system prompt, tools, or scorer made *after* looking at eval results is logged here with the reason. The goal is to make the agent behave better, not to make the numbers look better, so scorer changes are logged too.

## 2026-10-21: v1 baseline

- System prompt: `agent/system_prompt.md` as committed in Phase 4 (unchanged since).
- Model: `qwen3.5:9b` via local Ollama, `think=true`, temperature 0, seed 42.
- 15 questions (`eval/questions.yaml`), one fresh session per question.

**Issues seen in Phase 4 demos, before this eval** (not fixed, so the baseline shows them):
1. A top-5 list was reordered relative to the tool's padj ranking (GPX3/MAOA swap).
2. A QC-blocked answer did not cite the QC run.
3. A QC-blocked answer said "Please confirm if you want to proceed", which implies chat can override the gate (it can't: only the CLI flag or UI checkbox does).

**Scorer changes made while building the harness** (checked against the Phase 4 demo transcripts, before any live eval run):
- Required-number checks accept spelled-out numbers ("Thirteen"). The claim extractor still only reads digits; this blind spot is noted in the report.
- List items with no citation of their own inherit the list's citation (from a cited lead-in line ending in ":" or a citation-only line after the list). The report shows both this coverage and the strict per-line coverage.
- The sentence splitter no longer treats list markers ("1.") as sentences.

**v1 result** (`eval/report_v1_baseline.md`, transcripts `eval/transcripts/run_20261021T154517Z/`): 7/15 passed, 0/76 unsupported claims, 0 fabricated citations, agreement 2/7, refusal checks 3/5.

## 2026-10-21: v2, two fixes after reviewing v1 failures

### Agent fix: empty answers (affected 7/15 v1 answers: q01–q04, q06, q13, q14)
- **What happened:** with `think=true`, qwen3.5 sometimes wrote its final answer only into the hidden *thinking* channel and returned empty visible content. For q01 the hidden text was a correct, cited FKBP5 answer, but a user would have seen a blank reply.
- **Change (`agent/agent.py`):** if the model returns no tool calls *and* empty content, the agent sends one follow-up message ("Your last reply was empty. Write your final answer to the user now, following the rules and citing run ids.") and continues. A second empty reply gives an honest "I wasn't able to produce an answer" message. The retry is recorded in the transcript as a `retry` event.
- **Deliberately not done:** showing the thinking text as the answer. Thinking is not written for the user and isn't held to the citation rules.
- **System prompt:** unchanged.

### Scorer fix: comparison check (q05)
- **What happened:** q05 answered correctly ("KLF15 changes more … log2 fold change of 4.39 compared to PER1"), but the check only accepted "than", and its `[^.]*` pattern stopped at the decimal point in "4.39".
- **Change (`eval/run_eval.py`):** the check now works sentence by sentence and accepts "than / compared to / compared with / vs / versus / relative to". Tests added in `tests/test_hallucination.py`, including negative cases.
- **Effect isolated:** re-scoring the *unchanged* v1 transcripts with the fixed scorer gives 8/15 (`eval/report_v1_rescored.md`). Any improvement beyond that in v2 comes from the agent fix.

**v2 result** (`eval/report.md`, transcripts `eval/transcripts/run_20261021T204031Z/`): 14/15 passed, 2/120 unsupported claims (1 answer), 0 fabricated citations, agreement 7/7, refusal checks 5/5, citation support 29/29.
The empty-reply retry fired on exactly the 7 questions that were empty in v1 (q01–q04, q06, q13, q14), and all 7 then answered correctly. At temperature 0 with a fixed seed the empty-reply behavior was reproducible question by question.

**Open issues, not fixed (system prompt still unchanged):**
1. q13 suggests "BRCA1 or BRCA2" for the nonexistent BRCA7 without looking them up. Under the "never name a gene not in tool output" rule this counts as 2 unsupported claims.
2. Some restating sentences go uncited (q09 top-10 list, q11, q14). Every cited sentence is supported, but strict coverage is 56%.
3. The QC-blocked answer (q15) still doesn't cite the QC run and still says "Please confirm if you want to proceed". The scorer doesn't catch either one: the answer has no gene or number for the extractor to see. That's a scorer blind spot as well as an agent issue.
4. The ranking error seen in Phase 4 (top-5 *up* list, GPX3/MAOA swap) didn't recur on the top-5 *down* question in v1 or v2. With one sample per question, that's not evidence it's fixed.
