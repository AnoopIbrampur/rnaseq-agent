"""Tool-use agent: a local LLM (Ollama) that can only answer through pipeline tools.

The model never sees raw data. It sees compact tool results, each tagged with
the run_id of the provenance record that produced it. Three safeguards are
enforced in code rather than trusted to the prompt:

* QC gate: expression tools refuse to run until run_qc has passed in this
  session, unless the *human* set an override (logged as a run record).
* Citation check: every [run:...] in the answer must be a run from this
  session. Unknown ids are flagged as fabricated citations.
* Transcript: every question, tool call, tool result, and answer is saved,
  so the evaluation can re-check answers without calling the model again.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from agent.tools import DataSource, Session
from pipeline import config, provenance

SYSTEM_PROMPT = (Path(__file__).parent / "system_prompt.md").read_text()
TRANSCRIPTS = config.ROOT / "transcripts"
CITATION_RE = re.compile(r"\[run:([0-9a-f]{8})\]")
MAX_STEPS = 8
EMPTY_REPLY_NUDGE = ("Your last reply was empty. Write your final answer to the user now, "
                     "following the rules and citing run ids.")

TOOL_SCHEMAS: list[dict[str, Any]] = [
    {"type": "function", "function": {
        "name": "run_qc",
        "description": "Quality control of the 8 samples. Must be called first. Returns "
                       "passed (true/false), any flags, and a summary.",
        "parameters": {"type": "object", "properties": {}, "required": []}}},
    {"type": "function", "function": {
        "name": "run_de",
        "description": "Differential expression, dexamethasone vs untreated (design ~cell + dex). "
                       "Returns counts of significantly up/down genes and the top 10 genes "
                       "each way with log2 fold change and padj.",
        "parameters": {"type": "object", "properties": {}, "required": []}}},
    {"type": "function", "function": {
        "name": "run_enrichment",
        "description": "Pathway analysis (prerank GSEA, MSigDB Hallmark). Returns the gene sets "
                       "enriched among dexamethasone-induced and -repressed genes, with NES, "
                       "FDR, and leading-edge genes.",
        "parameters": {"type": "object", "properties": {}, "required": []}}},
    {"type": "function", "function": {
        "name": "lookup_gene",
        "description": "One gene's differential expression result and gene-set membership, or "
                       "the reason it was not tested. Accepts an official symbol or alias.",
        "parameters": {"type": "object",
                       "properties": {"symbol": {"type": "string",
                                                 "description": "Gene symbol, e.g. FKBP5"}},
                       "required": ["symbol"]}}},
]


# --- compact tool results for the model ----------------------------------------


def _fmt(x: Any) -> Any:
    """Round floats for the model: 4 significant digits, scientific when tiny."""
    if isinstance(x, float):
        return float(f"{x:.3g}") if abs(x) < 1e-3 and x != 0 else round(x, 4)
    if isinstance(x, dict):
        return {k: _fmt(v) for k, v in x.items()}
    if isinstance(x, list):
        return [_fmt(v) for v in x]
    return x


def compact_result(tool: str, rec: dict[str, Any]) -> dict[str, Any]:
    """What the model sees: the run_id plus the essential fields of the result."""
    r = rec["result"]
    out: dict[str, Any] = {"run_id": rec["run_id"], "summary": r.get("summary")}
    if tool == "run_qc":
        out |= {"passed": r["passed"], "flags": [f["message"] for f in r["flags"]]}
    elif tool == "run_de":
        keep = ["symbol", "log2FC", "padj"]
        out |= {"n_sig_up": r["n_sig_up"], "n_sig_down": r["n_sig_down"],
                "thresholds": "padj < 0.05 and |log2FC| > 1 (log2FC is shrunken)",
                "top_up": [{k: g[k] for k in keep} for g in r["top_up"][:10]],
                "top_down": [{k: g[k] for k in keep} for g in r["top_down"][:10]],
                "warnings": r["warnings"]}
    elif tool == "run_enrichment":
        def slim(p: dict[str, Any]) -> dict[str, Any]:
            return {"term": p["term"], "nes": p["nes"], "fdr": p["fdr"],
                    "leading_edge_first8": p["leading_edge"][:8]}
        out |= {"n_sig_up": r["n_sig_up"], "n_sig_down": r["n_sig_down"],
                "fdr_threshold": 0.05,
                "enriched_in_dex_induced": [slim(p) for p in r["top_up"]],
                "enriched_in_dex_repressed": [slim(p) for p in r["top_down"]]}
    elif tool == "lookup_gene":
        out |= {k: v for k, v in r.items() if k not in ("artifacts", "summary", "query")}
        out["gene_sets"] = [s for s in r.get("gene_sets", []) if s["significant"]]
    return _fmt(out)


# --- backends --------------------------------------------------------------------


class OllamaBackend:
    """Local model through the Ollama server."""

    def __init__(self, model: str, host: str | None = None, think: bool = True) -> None:
        import ollama
        self.client = ollama.Client(host=host) if host else ollama.Client()
        self.model, self.think = model, think

    def chat(self, messages: list[dict[str, Any]]) -> dict[str, Any]:
        try:
            r = self.client.chat(model=self.model, messages=messages, tools=TOOL_SCHEMAS,
                                 think=self.think,
                                 options={"num_ctx": 16384, "temperature": 0,
                                          "seed": config.SEED})
        except Exception as e:  # connection refused, model missing, ...
            raise RuntimeError(
                f"Ollama call failed ({e}). Is the server running (`ollama serve`) and is "
                f"the model pulled (`ollama pull {self.model}`)?"
            ) from e
        calls = [{"name": c.function.name, "arguments": dict(c.function.arguments or {})}
                 for c in (r.message.tool_calls or [])]
        return {"content": r.message.content or "", "tool_calls": calls,
                "thinking": getattr(r.message, "thinking", None)}


def make_backend() -> OllamaBackend:
    """Backend from environment (.env): AGENT_BACKEND, AGENT_MODEL, OLLAMA_HOST, AGENT_THINK."""
    from dotenv import load_dotenv
    load_dotenv(config.ROOT / ".env")
    backend = os.getenv("AGENT_BACKEND", "ollama")
    if backend != "ollama":
        raise RuntimeError(
            f"AGENT_BACKEND={backend!r} is not implemented; this project runs on local "
            "Ollama. Set AGENT_BACKEND=ollama."
        )
    return OllamaBackend(model=os.getenv("AGENT_MODEL", "qwen3.5:9b"),
                         host=os.getenv("OLLAMA_HOST") or None,
                         think=os.getenv("AGENT_THINK", "true").lower() == "true")


# --- agent -------------------------------------------------------------------------


@dataclass
class Answer:
    """The agent's reply plus everything needed to audit it."""

    question: str
    text: str
    citations: list[str]
    unknown_citations: list[str]
    tool_calls: list[dict[str, Any]]
    steps: int
    duration_s: float
    transcript_path: str | None = None


@dataclass
class Agent:
    """A conversation: one Session (data + runs), one backend, a message history."""

    backend: Any = field(default_factory=make_backend)
    session: Session = field(default_factory=Session)
    override_qc: bool = False
    override_reason: str = ""
    messages: list[dict[str, Any]] = field(default_factory=list)
    log: list[dict[str, Any]] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.messages = [{"role": "system", "content": SYSTEM_PROMPT}]
        self._qc: dict[str, Any] | None = None
        self._override_logged = False
        self._transcript_path: Path | None = None

    # -- tool dispatch with the QC gate --------------------------------------------

    def _gate(self) -> str | None:
        """Return an error message if expression tools are not allowed yet."""
        if self._qc is None:
            return "QC has not been run in this session. Call run_qc first."
        if not self._qc["result"]["passed"]:
            if self.override_qc:
                if not self._override_logged:
                    provenance.log_event("qc_override", {
                        "qc_run_id": self._qc["run_id"],
                        "flags": [f["message"] for f in self._qc["result"]["flags"]],
                        "reason": self.override_reason or "(no reason given)",
                        "set_by": "user",
                    })
                    self._override_logged = True
                return None
            return ("QC FAILED, so interpretation is blocked. Explain the QC flags to the "
                    "user. Only the user can override this, outside the chat.")
        return None

    def call_tool(self, name: str, args: dict[str, Any]) -> dict[str, Any]:
        """Run one tool call; return the compact result (or an error) for the model."""
        try:
            if name == "run_qc":
                self._qc = self.session.run_qc()
                return compact_result(name, self._qc)
            if name not in {"run_de", "run_enrichment", "lookup_gene"}:
                return {"error": f"Unknown tool {name!r}."}
            if (err := self._gate()) is not None:
                return {"error": err}
            if name == "lookup_gene":
                symbol = str(args.get("symbol", "")).strip()
                if not symbol:
                    return {"error": "lookup_gene needs a 'symbol' argument."}
                rec = self.session.lookup_gene(symbol)
            else:
                rec = getattr(self.session, name)()
            return compact_result(name, rec)
        except Exception as e:  # surface tool failures to the model and the log
            return {"error": f"{type(e).__name__}: {e}"}

    # -- main loop ------------------------------------------------------------------

    def ask(self, question: str, save: bool = True) -> Answer:
        """Answer one question, calling tools as needed; validate citations."""
        t0 = time.perf_counter()
        self.messages.append({"role": "user", "content": question})
        turn: dict[str, Any] = {"question": question, "events": []}
        calls: list[dict[str, Any]] = []
        text, steps, retried = "", 0, False
        for steps in range(1, MAX_STEPS + 1):
            reply = self.backend.chat(self.messages)
            turn["events"].append({"type": "model", "content": reply["content"],
                                   "tool_calls": reply["tool_calls"],
                                   "thinking": reply.get("thinking")})
            if not reply["tool_calls"] and not reply["content"].strip():
                # With thinking on, the model sometimes writes its answer only in the
                # hidden thinking channel. We never show thinking text to the user;
                # instead we ask once for the answer, then give up honestly.
                if retried:
                    text = ("I wasn't able to produce an answer for this question. "
                            "Please try rephrasing it.")
                    break
                retried = True
                self.messages.append({"role": "assistant", "content": ""})
                self.messages.append({"role": "user", "content": EMPTY_REPLY_NUDGE})
                turn["events"].append({"type": "retry", "reason": "empty_content",
                                       "message": EMPTY_REPLY_NUDGE})
                continue
            if not reply["tool_calls"]:
                text = reply["content"].strip()
                self.messages.append({"role": "assistant", "content": text})
                break
            self.messages.append({"role": "assistant", "content": reply["content"],
                                  "tool_calls": [{"function": c} for c in reply["tool_calls"]]})
            for c in reply["tool_calls"]:
                result = self.call_tool(c["name"], c["arguments"])
                calls.append({"name": c["name"], "arguments": c["arguments"],
                              "run_id": result.get("run_id"), "error": result.get("error")})
                turn["events"].append({"type": "tool", "name": c["name"],
                                       "arguments": c["arguments"], "result": result})
                self.messages.append({"role": "tool", "tool_name": c["name"],
                                      "content": json.dumps(result)})
        else:
            text = ("I could not finish within the tool-call limit, so I won't give an "
                    "answer I can't back up.")

        cited = list(dict.fromkeys(CITATION_RE.findall(text)))
        unknown = [c for c in cited if c not in self.session.runs]
        ans = Answer(question, text, cited, unknown, calls, steps,
                     round(time.perf_counter() - t0, 2))
        turn |= {"answer": text, "citations": cited, "unknown_citations": unknown,
                 "duration_s": ans.duration_s}
        self.log.append(turn)
        if save:
            ans.transcript_path = str(self.save_transcript())
        return ans

    def save_transcript(self, path: Path | None = None) -> Path:
        """Write the conversation, every tool call/result, and cited run records."""
        TRANSCRIPTS.mkdir(exist_ok=True)
        if path is None:
            if self._transcript_path is None:
                stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
                self._transcript_path = TRANSCRIPTS / f"{stamp}_{os.urandom(2).hex()}.json"
            path = self._transcript_path
        doc = {
            "backend": {"type": type(self.backend).__name__,
                        "model": getattr(self.backend, "model", None),
                        "think": getattr(self.backend, "think", None)},
            "data_source": {"counts": provenance.rel(self.session.source.counts),
                            "metadata": provenance.rel(self.session.source.metadata)},
            "override_qc": self.override_qc,
            "system_prompt_sha256": hashlib.sha256(SYSTEM_PROMPT.encode()).hexdigest(),
            "turns": self.log,
            "runs": {rid: {k: v for k, v in rec.items() if k != "result"} | {
                "result": rec["result"]} for rid, rec in self.session.runs.items()},
        }
        path.write_text(json.dumps(doc, indent=2, default=str))
        return path


def new_agent(counts: Path | None = None, metadata: Path | None = None,
              override_qc: bool = False, override_reason: str = "") -> Agent:
    """Convenience constructor used by the CLI and UI."""
    src = DataSource(counts or config.COUNTS_FILE, metadata or config.SERIES_MATRIX_FILE)
    return Agent(session=Session(source=src), override_qc=override_qc,
                 override_reason=override_reason)
