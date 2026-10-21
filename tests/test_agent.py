"""Agent loop tests with a scripted fake model (no Ollama needed)."""

from __future__ import annotations

import json

import pytest

from agent import agent as agent_mod
from agent.agent import Agent, compact_result
from agent.tools import DataSource, Session
from pipeline import config, data


class ScriptedBackend:
    """Replays a fixed list of replies; the 'answer' may reference run ids via {tool}."""

    model, think = "scripted", False

    def __init__(self, steps: list[dict]) -> None:
        self.steps, self.seen = steps, []

    def chat(self, messages: list[dict]) -> dict:
        self.seen.append(messages[-1])
        step = self.steps.pop(0)
        if "answer" in step:
            ids = {}
            for m in messages:
                if m["role"] == "tool":
                    r = json.loads(m["content"])
                    if "run_id" in r:
                        ids[m["tool_name"]] = r["run_id"]
            return {"content": step["answer"].format(**ids), "tool_calls": []}
        return {"content": "", "tool_calls": step["calls"]}


@pytest.fixture(autouse=True)
def tmp_transcripts(tmp_path, monkeypatch):
    monkeypatch.setattr(agent_mod, "TRANSCRIPTS", tmp_path / "transcripts")


def call(name: str, **args) -> dict:
    return {"calls": [{"name": name, "arguments": args}]}


def make(steps: list[dict], **kw) -> Agent:
    return Agent(backend=ScriptedBackend(steps), session=kw.pop("session", Session()), **kw)


def test_full_loop_with_valid_citation() -> None:
    a = make([call("run_qc"), call("lookup_gene", symbol="FKBP5"),
              {"answer": "FKBP5 is up [run:{lookup_gene}]."}])
    ans = a.ask("Is FKBP5 up?")
    assert ans.unknown_citations == [] and len(ans.citations) == 1
    assert [c["name"] for c in ans.tool_calls] == ["run_qc", "lookup_gene"]
    doc = json.loads(open(ans.transcript_path).read())
    assert doc["turns"][0]["answer"] == ans.text
    assert ans.citations[0] in doc["runs"]


def test_fabricated_citation_is_flagged() -> None:
    a = make([call("run_qc"), {"answer": "Everything is up [run:deadbeef]."}])
    ans = a.ask("?")
    assert ans.unknown_citations == ["deadbeef"]


def test_expression_tools_blocked_before_qc() -> None:
    a = make([call("run_de"), {"answer": "ok"}])
    a.ask("?")
    tool_event = a.log[0]["events"][1]
    assert "Call run_qc first" in tool_event["result"]["error"]


def test_failed_qc_blocks_and_override_is_logged(tmp_path, monkeypatch) -> None:
    bad = data.write_corrupted_counts("low_library")
    src = DataSource(counts=bad)
    blocked = make([call("run_qc"), call("run_de"), {"answer": "blocked"}],
                   session=Session(source=src))
    blocked.ask("?")
    assert "QC FAILED" in blocked.log[0]["events"][3]["result"]["error"]

    logged = []
    monkeypatch.setattr(agent_mod.provenance, "log_event",
                        lambda kind, details: logged.append((kind, details)))
    over = make([call("run_qc"), call("run_de"), {"answer": "done"}],
                session=Session(source=src), override_qc=True, override_reason="testing")
    over.ask("?")
    de_event = over.log[0]["events"][3]
    assert "error" not in de_event["result"] and de_event["result"]["n_sig_up"] >= 0
    assert logged and logged[0][0] == "qc_override"
    assert logged[0][1]["reason"] == "testing" and logged[0][1]["set_by"] == "user"


def test_tool_errors_are_returned_not_raised() -> None:
    a = make([call("run_qc"), call("lookup_gene"), call("bogus_tool"), {"answer": "x"}])
    a.ask("?")
    errors = [e["result"].get("error") for e in a.log[0]["events"] if e["type"] == "tool"]
    assert "needs a 'symbol'" in errors[1] and "Unknown tool" in errors[2]


def test_step_limit_refuses_to_answer(monkeypatch) -> None:
    monkeypatch.setattr(agent_mod, "MAX_STEPS", 2)
    a = make([call("run_qc"), call("run_qc")])
    ans = a.ask("?")
    assert "won't give an answer" in ans.text


def test_compact_result_keeps_run_id_and_rounds() -> None:
    rec = Session().lookup_gene("FKBP5")
    out = compact_result("lookup_gene", rec)
    assert out["run_id"] == rec["run_id"]
    assert out["symbol"] == "FKBP5" and out["direction"] == "up"
    assert all(s["significant"] for s in out["gene_sets"])


def test_system_prompt_example_uses_no_real_gene() -> None:
    assert "GENEX" in agent_mod.SYSTEM_PROMPT and "FKBP5" not in agent_mod.SYSTEM_PROMPT
    assert config.SEED == 42


def test_empty_reply_gets_one_retry_then_answer() -> None:
    a = make([call("run_qc"), {"answer": ""}, {"answer": "QC passed [run:{run_qc}]."}])
    ans = a.ask("Did QC pass?")
    assert ans.text.startswith("QC passed") and ans.unknown_citations == []
    assert any(e["type"] == "retry" for e in a.log[0]["events"])


def test_two_empty_replies_give_honest_fallback() -> None:
    a = make([call("run_qc"), {"answer": ""}, {"answer": "  "}])
    ans = a.ask("?")
    assert "wasn't able to produce an answer" in ans.text
