#!/usr/bin/env bash
# Terminal demo that needs no LLM: a recorded cited answer, its run record,
# a live reproduce of the DE run, and the eval headline table.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=.venv/bin/python
pause() { read -r -p $'\n[enter to continue] ' _; }

echo "== 1. Run the pipeline (cached after the first time) =="
$PY cli.py analyze --top 5 2>/dev/null | grep -v "WARNING\|order of those"
pause

echo "== 2. A recorded agent answer (eval transcript q01, local qwen3.5:9b) =="
T=$(ls -d eval/transcripts/run_* | tail -1)/q01.json
$PY - "$T" <<'PY'
import json, sys
d = json.load(open(sys.argv[1])); t = d["turns"][-1]
print("Q:", t["question"]); print("A:", t["answer"])
print("tools:", ", ".join(e["name"] for e in t["events"] if e["type"] == "tool"))
PY
pause

echo "== 3. The run record behind a DE claim =="
RID=$(grep -l '"tool": "run_de"' runs/*.json | xargs ls -t | head -1 | xargs basename | cut -d. -f1)
$PY - "$RID" <<'PY'
import json, sys
r = json.load(open(f"runs/{sys.argv[1]}.json")); r.pop("result")
print(json.dumps(r, indent=2)[:2500])
PY
pause

echo "== 4. Reproduce it from scratch =="
$PY cli.py reproduce "$RID" 2>/dev/null | grep -v "WARNING\|order of those"
pause

echo "== 5. Evaluation headline =="
sed -n '/## Headline numbers/,/## Per question/p' eval/report.md | sed '$d'
