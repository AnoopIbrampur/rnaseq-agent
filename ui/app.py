"""Streamlit chat UI: ask questions, inspect the provenance behind every claim.

Run:  .venv/bin/streamlit run ui/app.py

Design follows Apple's conventions: system fonts with size-specific tracking,
Messages-style bubbles for the user, grouped inset lists for records, a
translucent input bar, instant press feedback, and live status while the
local model works. Theme colors live in .streamlit/config.toml.
"""

from __future__ import annotations

import html
import json
import re
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from agent.agent import CITATION_RE, new_agent  # noqa: E402
from pipeline import data, provenance  # noqa: E402

st.set_page_config(page_title="Traceable RNA-seq Agent", page_icon="🧬", layout="centered",
                   initial_sidebar_state="expanded")

DATA_SOURCES = {
    "Airway study (real data)": None,
    "QC demo: scrambled sample": "scrambled_sample",
    "QC demo: under-sequenced sample": "low_library",
}
EXAMPLES = [
    "Which pathways change when these cells are treated with dexamethasone?",
    "Is FKBP5 affected, and how confident are we?",
    "How many genes are upregulated, and what are the top 5?",
]
TOOL_NAMES = {
    "run_qc": "Quality control",
    "run_de": "Differential expression",
    "run_enrichment": "Pathway enrichment",
    "lookup_gene": "Gene lookup",
    "qc_override": "QC override",
}
TOOL_PROGRESS = {
    "run_qc": "Checking data quality…",
    "run_de": "Running differential expression…",
    "run_enrichment": "Running pathway enrichment…",
    "lookup_gene": "Looking up {symbol}…",
}

# --- theme -------------------------------------------------------------------------

LIGHT = {
    "accent": "#007aff", "accent-fill": "rgba(0,122,255,0.12)", "on-accent": "#ffffff",
    "fill": "#f2f2f7", "cell": "#ffffff", "label": "#1d1d1f", "label-2": "rgba(60,60,67,0.62)",
    "sep": "rgba(60,60,67,0.18)", "green": "#248a3d", "green-fill": "rgba(52,199,89,0.14)",
    "orange": "#c93400", "orange-fill": "rgba(255,149,0,0.14)", "bar": "rgba(255,255,255,0.72)",
}
DARK = {
    "accent": "#0a84ff", "accent-fill": "rgba(10,132,255,0.18)", "on-accent": "#ffffff",
    "fill": "#1c1c1e", "cell": "#2c2c2e", "label": "#f5f5f7", "label-2": "rgba(235,235,245,0.62)",
    "sep": "rgba(84,84,88,0.65)", "green": "#30d158", "green-fill": "rgba(48,209,88,0.16)",
    "orange": "#ff9f0a", "orange-fill": "rgba(255,159,10,0.16)", "bar": "rgba(0,0,0,0.66)",
}


def css_vars(tokens: dict[str, str]) -> str:
    return "".join(f"--{k}:{v};" for k, v in tokens.items())


def theme_css() -> str:
    """Token block following the OS light/dark setting, as the Streamlit theme does."""
    return (f":root{{{css_vars(LIGHT)}}}"
            f"@media (prefers-color-scheme: dark){{:root{{{css_vars(DARK)}}}}}")


STYLE = """
/* chrome: translucent header, no decoration bar */
[data-testid="stDecoration"] { display: none; }
[data-testid="stHeader"] {
  background: var(--bar);
  -webkit-backdrop-filter: blur(20px) saturate(180%); backdrop-filter: blur(20px) saturate(180%);
}
[data-testid="stMainBlockContainer"] { padding-top: 3.5rem; padding-bottom: 7rem; }

/* type: size-specific tracking */
html, body, [data-testid="stAppViewContainer"] { -webkit-font-smoothing: antialiased; }
.stMarkdown p, .stMarkdown li { font-size: 1.0625rem; line-height: 1.47; letter-spacing: -0.01em; }
.hero h1 {
  font-size: clamp(2.1rem, 5.5vw, 2.6rem) !important; line-height: 1.08 !important;
  letter-spacing: -0.025em; font-weight: 700 !important; margin: 0 0 .35rem !important; padding: 0 !important;
}
.hero p { color: var(--label-2); font-size: 1.0625rem; letter-spacing: -0.01em; margin: 0 0 .9rem; max-width: 36rem; }
.section-label {
  font-size: .8125rem; letter-spacing: .01em; text-transform: uppercase;
  color: var(--label-2); margin: 1.5rem 0 .5rem; padding-left: 1rem;
}

/* status pills */
.pills { display: flex; flex-wrap: wrap; gap: .4rem; margin-bottom: .5rem; }
.pill {
  display: inline-flex; align-items: center; gap: .4rem;
  font-size: .8125rem; font-weight: 590; letter-spacing: -0.005em;
  padding: .3rem .7rem; border-radius: 999px; background: var(--fill); color: var(--label);
}
.pill .dot { width: .5rem; height: .5rem; border-radius: 50%; background: var(--label-2); }
.pill.ok { background: var(--green-fill); color: var(--green); } .pill.ok .dot { background: var(--green); }
.pill.warn { background: var(--orange-fill); color: var(--orange); } .pill.warn .dot { background: var(--orange); }

/* buttons: instant press feedback */
.stButton button { transition: transform 100ms ease-out, background-color 100ms ease-out; }
.stButton button:active { transform: scale(0.97); }

/* suggestions: an inset grouped list of tappable rows */
[class*="st-key-suggest"] button {
  width: 100%; justify-content: flex-start; text-align: left;
  border-radius: 0 !important; background: var(--fill); border: none; box-shadow: none;
  padding: .8rem 2.25rem .8rem 1rem; min-height: 0; position: relative; color: var(--label);
}
[class*="st-key-suggest"] button p { font-size: 1rem; letter-spacing: -0.01em; font-weight: 400; text-align: left; }
[class*="st-key-suggest"] button > div { justify-content: flex-start; width: 100%; }
[class*="st-key-suggest"] button::after {
  content: "›"; position: absolute; right: 1rem; top: 50%; transform: translateY(-52%);
  font-size: 1.35rem; color: var(--label-2);
}
[class*="st-key-suggest"] button:hover { background: var(--accent-fill); color: var(--label); }
[class*="st-key-suggest"] button:active { transform: none; background: var(--accent-fill); }
.st-key-suggest-0 button { border-radius: 12px 12px 0 0 !important; }
.st-key-suggest-2 button { border-radius: 0 0 12px 12px !important; }
.st-key-suggestions, .st-key-suggestions [data-testid="stVerticalBlock"] { gap: 0 !important; }
.st-key-suggest-1 button, .st-key-suggest-2 button { box-shadow: inset 1rem 0.5px 0 0 var(--fill), inset 0 0.5px 0 0 var(--sep); }

/* chat: Messages-style user bubble, open assistant text */
[data-testid="stChatMessage"] { background: transparent !important; padding: .25rem 0 !important; gap: 0; }
[data-testid="stChatMessageAvatarUser"], [data-testid="stChatMessageAvatarAssistant"] { display: none; }
[data-testid="stChatMessage"]:has([data-testid="stChatMessageAvatarUser"]) { justify-content: flex-end; }
[data-testid="stChatMessage"]:has([data-testid="stChatMessageAvatarUser"]) [data-testid="stChatMessageContent"] {
  flex: 0 1 auto; max-width: 80%; background: var(--accent); border-radius: 1.25rem;
  padding: .55rem .95rem; margin-left: auto !important; margin-right: 0 !important;
}
[data-testid="stChatMessage"]:has([data-testid="stChatMessageAvatarUser"]) [data-testid="stChatMessageContent"] * {
  color: var(--on-accent) !important;
}
.meta { font-size: .8125rem; color: var(--label-2); letter-spacing: 0; margin: .35rem 0 .6rem; }
.cite {
  font: 500 .75rem/1 ui-monospace, "SF Mono", Menlo, monospace; letter-spacing: 0;
  color: var(--accent); background: var(--accent-fill);
  padding: .18rem .4rem; border-radius: 6px; white-space: nowrap; vertical-align: .08em;
}
.cite.bad { color: var(--orange); background: var(--orange-fill); }
.notice {
  font-size: .9375rem; padding: .7rem 1rem; border-radius: 12px; margin: .4rem 0;
  background: var(--orange-fill); color: var(--orange);
}

/* provenance cards: expanders as inset grouped cards */
[data-testid="stExpander"] details {
  border: none !important; background: var(--fill); border-radius: 12px !important; overflow: hidden;
}
[data-testid="stExpander"] summary { padding: .7rem 1rem; }
[data-testid="stExpander"] summary:hover { background: var(--accent-fill); }
[data-testid="stExpander"] summary p { font-size: .9375rem; letter-spacing: -0.01em; }

/* grouped key/value lists */
.group { background: var(--cell); border-radius: 12px; overflow: hidden; margin: .25rem 0 .75rem; }
.kv { display: grid; grid-template-columns: minmax(6.5rem, 34%) 1fr; gap: .75rem; padding: .6rem 1rem; position: relative; }
.kv + .kv::before { content: ""; position: absolute; top: 0; left: 1rem; right: 0; height: .5px; background: var(--sep); }
.kv .k { color: var(--label-2); font-size: .875rem; }
.kv .v { font-size: .875rem; color: var(--label); min-width: 0; overflow-wrap: anywhere; }
.mono { font-family: ui-monospace, "SF Mono", Menlo, monospace; font-size: .8125rem; letter-spacing: 0; }
.summary-text { font-size: .9375rem; color: var(--label); line-height: 1.45; padding: .2rem .1rem .6rem; }

/* sidebar */
[data-testid="stSidebar"] h2 { font-size: 1.375rem !important; letter-spacing: -0.02em; }
.runs { background: var(--cell); border-radius: 12px; overflow: hidden; margin-bottom: 1rem; }
.runs .r { display: flex; justify-content: space-between; gap: .5rem; padding: .55rem .85rem; position: relative; font-size: .875rem; }
.runs .r + .r::before { content: ""; position: absolute; top: 0; left: .85rem; right: 0; height: .5px; background: var(--sep); }
.runs .t { color: var(--label); } .runs .s { color: var(--label-2); }
.runs .empty { padding: .7rem .85rem; color: var(--label-2); font-size: .875rem; }

/* input bar: translucent material, content scrolls beneath */
[data-testid="stBottom"] > div {
  background: var(--bar) !important;
  -webkit-backdrop-filter: blur(20px) saturate(180%); backdrop-filter: blur(20px) saturate(180%);
}
[data-testid="stChatInput"] { border-radius: 1.4rem !important; }

@media (prefers-reduced-motion: reduce) {
  .stButton button { transition: none; } .stButton button:active { transform: none; }
}
@media (prefers-reduced-transparency: reduce) {
  [data-testid="stHeader"], [data-testid="stBottom"] > div {
    -webkit-backdrop-filter: none; backdrop-filter: none; background: var(--fill) !important;
  }
}
@media (prefers-contrast: more) {
  .group, .runs, [data-testid="stExpander"] details, .pill { box-shadow: 0 0 0 1px var(--sep); }
}
"""


def inject_style() -> None:
    st.markdown(f"<style>{theme_css()}{STYLE}</style>", unsafe_allow_html=True)


# --- session ----------------------------------------------------------------------------


def counts_path(choice: str) -> Path | None:
    kind = DATA_SOURCES[choice]
    return data.write_corrupted_counts(kind) if kind else None


def reset_agent() -> None:
    s = st.session_state
    s.agent = new_agent(counts_path(s.source), override_qc=s.override,
                        override_reason=s.get("reason", ""))
    s.history = []


def qc_state() -> tuple[str, str]:
    """(pill class, label) for the latest QC run in this conversation."""
    runs = st.session_state.agent.session.runs
    qc = [r for r in runs.values() if r["tool"] == "run_qc"]
    if not qc:
        return "", "QC runs with your first question"
    if qc[-1]["result"]["passed"]:
        return "ok", "QC passed"
    return "warn", "QC failed · interpretation blocked" if not st.session_state.override \
        else "QC failed · override on"


# --- rendering helpers -------------------------------------------------------------------


def esc(x: object) -> str:
    return html.escape(str(x), quote=True)


def kv_group(rows: list[tuple[str, str]]) -> str:
    """Grouped key/value list; values are already-escaped HTML."""
    return '<div class="group">' + "".join(
        f'<div class="kv"><span class="k">{esc(k)}</span><span class="v">{v}</span></div>'
        for k, v in rows) + "</div>"


def answer_html(text: str, known: set[str]) -> str:
    """Escape model text, then turn [run:id] into citation pills."""
    def pill(m: re.Match) -> str:
        rid = m.group(1)
        if rid in known:
            return f'<span class="cite" title="Provenance record below">run:{rid}</span>'
        return f'<span class="cite bad" title="Not a run in this conversation">run:{rid} ⚠︎</span>'
    return CITATION_RE.sub(pill, html.escape(text, quote=False))


def mentioned_symbols(text: str, table: pd.DataFrame) -> list[str]:
    """Gene symbols from the DE table that appear as words in the answer."""
    words = set(re.findall(r"\b[A-Z][A-Z0-9-]{1,14}\b", text))
    return [s for s in table["symbol"].dropna().unique() if s in words]


def short(sha: str, n: int = 12) -> str:
    return f'<span class="mono">{esc(sha[:n])}…</span>'


def render_record(rec: dict, answer_text: str) -> None:
    """Provenance card for one run: summary, inputs/outputs, parameters, data."""
    tool = TOOL_NAMES.get(rec["tool"], rec["tool"])
    status = "Reused from cache" if rec.get("cache_hit") else "Computed"
    with st.expander(f"run:{rec['run_id']}  ·  {tool}  ·  {status}"):
        t_sum, t_io, t_par, t_data = st.tabs(["Summary", "Inputs & outputs", "Parameters", "Data"])
        outs = {k: provenance.absolute(v) for k, v in rec["outputs"].items()}
        with t_sum:
            st.markdown(f'<div class="summary-text">{esc(rec.get("summary") or "")}</div>',
                        unsafe_allow_html=True)
            git = rec["git"]
            rows = [
                ("Run", f'<span class="mono">{esc(rec["run_id"])}</span>'),
                ("Status", f'{esc(status)} in {rec.get("duration_s", 0)} s'),
                ("Recorded", esc(rec["timestamp_utc"].replace("T", " ").replace("+00:00", " UTC"))),
                ("Git commit", f'<span class="mono">{esc((git["commit"] or "")[:10])}</span>'
                               + (" · uncommitted changes" if git["pipeline_dirty"] else "")),
                ("Tool code", short(rec["code_sha256"])),
            ]
            if rec["upstream_runs"]:
                rows.append(("Built on", " ".join(f'<span class="cite">run:{esc(r)}</span>'
                                                  for r in rec["upstream_runs"])))
            st.markdown(kv_group(rows), unsafe_allow_html=True)
            st.caption("Re-run from scratch and verify the outputs match:")
            st.code(f".venv/bin/python cli.py reproduce {rec['run_id']}", language="bash")
        with t_io:
            st.markdown('<div class="section-label" style="margin-top:.5rem">Inputs</div>'
                        + kv_group([(k, f'{esc(v["path"])}<br>{short(v["sha256"], 16)}')
                                    for k, v in rec["inputs"].items()]),
                        unsafe_allow_html=True)
            st.markdown('<div class="section-label">Outputs</div>'
                        + kv_group([(k, f'{esc(v)}<br>{short(rec["output_sha256"][k], 16)}')
                                    for k, v in rec["outputs"].items()]),
                        unsafe_allow_html=True)
            sw = rec["software"]
            st.markdown('<div class="section-label">Software</div>'
                        + kv_group([(k, f'<span class="mono">{esc(v)}</span>') for k, v in sw.items()]),
                        unsafe_allow_html=True)
        with t_par:
            st.json(rec["parameters"], expanded=1)
        with t_data:
            if rec["tool"] == "run_de":
                t = pd.read_csv(outs["table"])
                syms = mentioned_symbols(answer_text, t)
                rows = t[t["symbol"].isin(syms)] if syms else t[t["significant"]].head(10)
                st.caption("Rows behind this answer" if syms else "Top significant genes")
                st.dataframe(rows[["symbol", "baseMean", "log2FC", "lfcSE", "stat", "padj",
                                   "direction"]], hide_index=True, width="stretch")
            elif rec["tool"] == "run_enrichment":
                t = pd.read_csv(outs["table"])
                st.dataframe(t[["term", "nes", "fdr", "set_size", "leading_edge"]].head(15),
                             hide_index=True, width="stretch")
            elif rec["tool"] == "lookup_gene":
                st.json(json.loads(outs["lookup"].read_text()), expanded=1)
            elif rec["tool"] == "run_qc":
                r = rec["result"]
                flags = [f["message"] for f in r["flags"]] or ["No flags"]
                st.markdown(kv_group([("Passed", "Yes" if r["passed"] else "No")]
                                     + [("Flag", esc(f)) for f in flags]), unsafe_allow_html=True)
            for k, p in outs.items():
                if p.suffix == ".png":
                    st.image(str(p), caption=k)


def render_answer(turn: dict, known: set[str]) -> None:
    st.markdown(answer_html(turn["answer"], known), unsafe_allow_html=True)
    if turn["unknown_citations"]:
        st.markdown(f'<div class="notice">Citation check failed: '
                    f'{esc(", ".join(turn["unknown_citations"]))} '
                    f'{"is" if len(turn["unknown_citations"]) == 1 else "are"} not runs from this '
                    f'conversation.</div>', unsafe_allow_html=True)
    tools = [e for e in turn["events"] if e["type"] == "tool"]
    trail = " → ".join(
        TOOL_NAMES.get(e["name"], e["name"]) + (f" ({esc(e['arguments']['symbol'])})"
                                               if e["arguments"].get("symbol") else "")
        for e in tools) or "no new tools (used earlier results)"
    st.markdown(f'<div class="meta">{turn["duration_s"]:.0f} s · {trail}</div>',
                unsafe_allow_html=True)
    runs = st.session_state.agent.session.runs
    for rid in turn["citations"]:
        if rid in runs:
            render_record(runs[rid], turn["answer"])
    uncited = [e["result"]["run_id"] for e in tools
               if e["result"].get("run_id") and e["result"]["run_id"] not in turn["citations"]]
    if uncited:
        st.markdown(f'<div class="meta">Also used, not cited: '
                    + " ".join(f'<span class="cite">run:{esc(r)}</span>' for r in uncited)
                    + "</div>", unsafe_allow_html=True)


def render_turn(turn: dict) -> None:
    with st.chat_message("user"):
        st.markdown(esc(turn["question"]))
    with st.chat_message("assistant"):
        render_answer(turn, set(st.session_state.agent.session.runs))


# --- page -------------------------------------------------------------------------------------

inject_style()
s = st.session_state
s.setdefault("source", list(DATA_SOURCES)[0])
s.setdefault("override", False)
if "agent" not in s:
    reset_agent()

with st.sidebar:
    st.markdown("## Conversation")
    st.selectbox("Dataset", list(DATA_SOURCES), key="source", on_change=reset_agent)
    st.toggle("Interpret even if QC fails", key="override", on_change=reset_agent,
              help="Only you can turn this on; the agent can't. Each override is saved to runs/ "
                   "as its own record.")
    if s.override:
        st.text_input("Reason for overriding QC", key="reason", on_change=reset_agent,
                      placeholder="e.g. checking the flagged sample by eye")
    if st.button("New conversation", width="stretch"):
        reset_agent()
        st.rerun()
    st.markdown('<div class="section-label" style="padding-left:.85rem">Runs in this conversation</div>',
                unsafe_allow_html=True)
    runs = s.agent.session.runs
    rows = "".join(
        f'<div class="r"><span class="t"><span class="mono">{esc(rid)}</span> · '
        f'{esc(TOOL_NAMES.get(r["tool"], r["tool"]))}</span>'
        f'<span class="s">{"cache" if r.get("cache_hit") else "new"}</span></div>'
        for rid, r in runs.items()) or '<div class="empty">None yet</div>'
    st.markdown(f'<div class="runs">{rows}</div>', unsafe_allow_html=True)
    st.caption(f"Model: {getattr(s.agent.backend, 'model', '?')}, running on this computer "
               "with Ollama.")

qc_cls, qc_label = qc_state()
st.markdown(
    '<div class="hero"><h1>Ask the airway study</h1>'
    "<p>Dexamethasone vs untreated airway smooth muscle cells, 4 donors (GSE52778). "
    "Answers come only from analyses run here, and every claim links to its run record.</p></div>"
    f'<div class="pills"><span class="pill">{esc(s.source)}</span>'
    f'<span class="pill {qc_cls}"><span class="dot"></span>{esc(qc_label)}</span></div>',
    unsafe_allow_html=True)

for turn in s.history:
    render_turn(turn)

# Read the input first so the suggestions disappear the moment a question is sent.
question = st.chat_input("Ask about the dexamethasone response…") or s.pop("pending", None)

if not s.history and not question:
    st.markdown('<div class="section-label">Try asking</div>', unsafe_allow_html=True)
    with st.container(key="suggestions"):
        for i, q in enumerate(EXAMPLES):
            if st.button(q, key=f"suggest-{i}", width="stretch"):
                s.pending = q
                st.rerun()

if question:
    with st.chat_message("user"):
        st.markdown(esc(question))
    with st.chat_message("assistant"):
        status = st.status("Reading your question…", expanded=False)

        def progress(event: dict) -> None:
            if event["stage"] == "tool_start":
                label = TOOL_PROGRESS.get(event["name"], "Working…")
                status.update(label=label.format(symbol=event["arguments"].get("symbol", "the gene")))
                status.write(f"{TOOL_NAMES.get(event['name'], event['name'])}…")
            elif event["stage"] == "tool_done":
                r = event["result"]
                status.write(f"✓ run:{r['run_id']}" if "run_id" in r else f"✕ {r.get('error', '')}")
            elif event["stage"] == "model" and event["step"] > 1:
                status.update(label="Writing the answer…")

        try:
            ans = s.agent.ask(question, on_progress=progress)
            s.history.append(s.agent.log[-1])
            status.update(label=f"Answered in {ans.duration_s:.0f} s", state="complete")
        except RuntimeError as e:
            status.update(label="Couldn't reach the local model", state="error")
            st.markdown(f'<div class="notice">{esc(e)}</div>', unsafe_allow_html=True)
            st.stop()
    st.rerun()
