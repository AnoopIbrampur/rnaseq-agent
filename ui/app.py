"""Streamlit chat UI: ask questions, inspect the provenance behind every claim.

Run:  .venv/bin/streamlit run ui/app.py
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from agent.agent import CITATION_RE, new_agent  # noqa: E402
from pipeline import data, provenance  # noqa: E402

st.set_page_config(page_title="Traceable RNA-seq Agent", page_icon="🧬", layout="wide")

DATA_SOURCES = {
    "Airway GSE52778 (real data)": None,
    "Corrupted copy: one scrambled sample (QC demo)": "scrambled_sample",
    "Corrupted copy: one under-sequenced sample (QC demo)": "low_library",
}
EXAMPLES = [
    "Which pathways change when these cells are treated with dexamethasone?",
    "Is FKBP5 affected, and how confident are we?",
    "How many genes are upregulated, and what are the top 5?",
]


def counts_path(choice: str) -> Path | None:
    kind = DATA_SOURCES[choice]
    return data.write_corrupted_counts(kind) if kind else None


def reset_agent() -> None:
    s = st.session_state
    s.agent = new_agent(counts_path(s.source), override_qc=s.override,
                        override_reason=s.reason)
    s.history = []


def render_answer(text: str, known: set[str]) -> str:
    """Turn [run:id] into inline code badges; mark unknown ids loudly."""
    def badge(m: re.Match) -> str:
        rid = m.group(1)
        return f"`run:{rid}`" if rid in known else f"**:red[run:{rid} — NOT A RUN IN THIS SESSION]**"
    return CITATION_RE.sub(badge, text)


def mentioned_symbols(text: str, table: pd.DataFrame) -> list[str]:
    """Gene symbols from the DE table that appear as words in the answer."""
    words = set(re.findall(r"\b[A-Z][A-Z0-9-]{1,14}\b", text))
    return [s for s in table["symbol"].dropna().unique() if s in words]


def render_record(rec: dict, answer_text: str) -> None:
    """Expandable provenance panel for one run record."""
    with st.expander(f"run:{rec['run_id']} · {rec['tool']} · "
                     f"{'cache hit' if rec.get('cache_hit') else 'computed'} · "
                     f"{rec.get('duration_s', 0)}s"):
        st.caption(rec.get("summary") or "")
        c1, c2 = st.columns(2)
        with c1:
            st.markdown("**Inputs**")
            st.dataframe(pd.DataFrame([{"name": k, "path": v["path"], "sha256": v["sha256"][:16] + "…"}
                                       for k, v in rec["inputs"].items()]),
                         hide_index=True, width="stretch")
            st.markdown("**Parameters**")
            st.json(rec["parameters"], expanded=False)
        with c2:
            st.markdown("**Software & code**")
            st.json({**rec["software"], "code_sha256": rec["code_sha256"][:16] + "…",
                     "git_commit": (rec["git"]["commit"] or "")[:10],
                     "pipeline_dirty": rec["git"]["pipeline_dirty"]}, expanded=False)
            st.markdown("**Outputs**")
            st.dataframe(pd.DataFrame([{"output": k, "path": v,
                                        "sha256": rec["output_sha256"][k][:16] + "…"}
                                       for k, v in rec["outputs"].items()]),
                         hide_index=True, width="stretch")
        if rec["upstream_runs"]:
            st.markdown(f"**Upstream runs:** {', '.join('run:' + r for r in rec['upstream_runs'])}")

        outs = {k: provenance.absolute(v) for k, v in rec["outputs"].items()}
        if rec["tool"] == "run_de":
            t = pd.read_csv(outs["table"])
            syms = mentioned_symbols(answer_text, t)
            rows = t[t["symbol"].isin(syms)] if syms else t[t["significant"]].head(10)
            st.markdown("**Rows behind this answer**" if syms else "**Top significant genes**")
            st.dataframe(rows[["symbol", "gene_id", "baseMean", "log2FC", "lfcSE", "stat",
                               "padj", "direction"]], hide_index=True, width="stretch")
        elif rec["tool"] == "run_enrichment":
            t = pd.read_csv(outs["table"])
            st.dataframe(t[["term", "nes", "fdr", "set_size", "leading_edge"]].head(15),
                         hide_index=True, width="stretch")
        elif rec["tool"] == "lookup_gene":
            st.json(json.loads(outs["lookup"].read_text()), expanded=True)
        elif rec["tool"] == "run_qc":
            st.json({"passed": rec["result"]["passed"],
                     "flags": [f["message"] for f in rec["result"]["flags"]]})
        for k, p in outs.items():
            if p.suffix == ".png":
                st.image(str(p), caption=k)
        st.code(f".venv/bin/python cli.py reproduce {rec['run_id']}", language="bash")


def render_turn(turn: dict) -> None:
    agent = st.session_state.agent
    with st.chat_message("user"):
        st.write(turn["question"])
    with st.chat_message("assistant"):
        st.markdown(render_answer(turn["answer"], set(agent.session.runs)))
        if turn["unknown_citations"]:
            st.error(f"Citation check failed: {', '.join(turn['unknown_citations'])} "
                     "are not runs from this session.")
        tools = [e for e in turn["events"] if e["type"] == "tool"]
        st.caption(f"{turn['duration_s']}s · tools: "
                   + ", ".join(f"{e['name']}({e['arguments'].get('symbol', '')})" for e in tools))
        for rid in turn["citations"]:
            if rid in agent.session.runs:
                render_record(agent.session.runs[rid], turn["answer"])
        uncited = [e["result"].get("run_id") for e in tools
                   if e["result"].get("run_id") and e["result"]["run_id"] not in turn["citations"]]
        if uncited:
            with st.expander(f"Other runs used but not cited ({len(uncited)})"):
                for rid in uncited:
                    st.markdown(f"- `run:{rid}` · {agent.session.runs[rid]['tool']}")


# --- page --------------------------------------------------------------------------

s = st.session_state
with st.sidebar:
    st.header("Session")
    st.selectbox("Data", list(DATA_SOURCES), key="source", on_change=reset_agent)
    st.checkbox("Override failed QC (logged)", key="override", on_change=reset_agent,
                help="Lets the agent interpret results even when QC fails. The override is "
                     "written to runs/ as its own record. Only you can set this.")
    st.text_input("Override reason", key="reason", on_change=reset_agent)
    if st.button("New session"):
        reset_agent()
    if "agent" in s:
        b = s.agent.backend
        st.caption(f"Model: {getattr(b, 'model', '?')} (local Ollama)")
        st.subheader("Runs this session")
        for rid, rec in s.agent.session.runs.items():
            st.caption(f"`run:{rid}` {rec['tool']}{' (cache)' if rec.get('cache_hit') else ''}")

if "agent" not in s:
    reset_agent()

st.title("🧬 Traceable RNA-seq interpretation agent")
st.caption("Dexamethasone vs untreated airway smooth muscle cells (GSE52778). The agent can only "
           "state what the analysis tools computed, and every claim links to its run record.")

for turn in s.history:
    render_turn(turn)

if not s.history:
    cols = st.columns(len(EXAMPLES))
    for col, q in zip(cols, EXAMPLES):
        if col.button(q, width="stretch"):
            s.pending = q
            st.rerun()

question = st.chat_input("Ask about this experiment…") or s.pop("pending", None)
if question:
    with st.spinner("Running tools and writing an answer (local model, ~20-60s)…"):
        try:
            s.agent.ask(question)
            s.history.append(s.agent.log[-1])
        except RuntimeError as e:
            st.error(str(e))
    st.rerun()
