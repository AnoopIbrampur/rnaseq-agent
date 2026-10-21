"""Provenance: every tool call produces a run record, and any record can be replayed.

`execute()` is the only way tools run. It:
1. checksums every input file,
2. builds a cache key from (tool, input checksums, parameters, package
   versions, tool source code) and returns the stored result on a hit,
3. otherwise computes into a content-addressed folder results/<tool>/<key12>/,
4. writes runs/<run_id>.json describing exactly what happened.

`reproduce()` re-executes a record from scratch (never from cache) in a temp
folder and checks that every output matches.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import subprocess
import tempfile
import time
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

from pipeline import cache, config, data

PACKAGES = ["pydeseq2", "gseapy", "pandas", "numpy", "scipy", "anndata", "matplotlib"]


# --- environment -------------------------------------------------------------


def software_versions() -> dict[str, str]:
    """Python and scientific package versions."""
    return {"python": platform.python_version(), **{p: version(p) for p in PACKAGES}}


def git_state() -> dict[str, Any]:
    """Current commit and whether the working tree has uncommitted changes."""
    def run(*args: str) -> str:
        return subprocess.run(["git", *args], cwd=config.ROOT, capture_output=True,
                              text=True, check=False).stdout.strip()
    commit = run("rev-parse", "HEAD") or None
    dirty = bool(run("status", "--porcelain", "--", "pipeline")) if commit else None
    return {"commit": commit, "pipeline_dirty": dirty}


def code_sha(modules: list[str]) -> str:
    """SHA-256 over the source of the pipeline modules a tool depends on."""
    h = hashlib.sha256()
    for m in sorted(modules):
        h.update((config.ROOT / "pipeline" / f"{m}.py").read_bytes())
    return h.hexdigest()


def rel(path: Path | str) -> str:
    """Path relative to the project root when possible (portable records)."""
    p = Path(path).resolve()
    try:
        return str(p.relative_to(config.ROOT))
    except ValueError:
        return str(p)


def absolute(path: str) -> Path:
    p = Path(path)
    return p if p.is_absolute() else config.ROOT / p


def to_jsonable(obj: Any) -> Any:
    """Convert numpy/pandas scalars and NaN so results round-trip through JSON."""
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, np.generic):
        obj = obj.item()
    if isinstance(obj, float) and math.isnan(obj):
        return None
    if isinstance(obj, Path):
        return rel(obj)
    return obj


# --- tool registry -------------------------------------------------------------


@dataclass(frozen=True)
class ToolSpec:
    """How to run a tool from file inputs, and which source files define it."""

    compute: Callable[[dict[str, Path], dict[str, Any], Path], dict[str, Any]]
    modules: list[str]


def _compute_qc(inputs: dict[str, Path], params: dict[str, Any], out: Path) -> dict[str, Any]:
    from pipeline import qc
    meta = data.load_metadata(inputs["metadata"])
    counts = data.load_counts(inputs["counts"], samples=list(meta.index))
    return qc.run_qc(counts, meta, out, rules=params["qc_rules"], low_count=params["low_count"])


def _compute_de(inputs: dict[str, Path], params: dict[str, Any], out: Path) -> dict[str, Any]:
    from pipeline import de
    ds = data.load_dataset(inputs["counts"], inputs["metadata"], inputs["gene_info"],
                           inputs["gene_history"])
    return de.run_de(ds.counts, ds.metadata, ds.gene_map, out,
                     params=params["de"], low_count=params["low_count"])


def _compute_enrichment(inputs: dict[str, Path], params: dict[str, Any], out: Path) -> dict[str, Any]:
    from pipeline import enrichment
    return enrichment.run_enrichment(pd.read_csv(inputs["de_table"]), inputs["gene_sets"],
                                     out, params=params)


def _compute_lookup(inputs: dict[str, Path], params: dict[str, Any], out: Path) -> dict[str, Any]:
    from pipeline import genes
    result = genes.lookup_gene(
        params["query"],
        de_table=pd.read_csv(inputs["de_table"]),
        filtered=pd.read_csv(inputs["filtered_genes"]),
        gsea_table=pd.read_csv(inputs["gsea_table"]),
        gene_sets=genes.read_gmt(inputs["gene_sets"]),
        gene_map=data.load_gene_map(inputs["gene_info"], inputs["gene_history"]),
        fdr_threshold=params["fdr_threshold"],
        alpha=params["alpha"],
        lfc_threshold=params["lfc_threshold"],
        min_samples=params["min_samples"],
    )
    out.mkdir(parents=True, exist_ok=True)
    path = out / "lookup.json"
    path.write_text(json.dumps(to_jsonable(result), indent=2))
    return {**result, "artifacts": {"lookup": str(path)}}


TOOLS: dict[str, ToolSpec] = {
    "run_qc": ToolSpec(_compute_qc, ["qc", "data"]),
    "run_de": ToolSpec(_compute_de, ["de", "qc", "data"]),
    "run_enrichment": ToolSpec(_compute_enrichment, ["enrichment"]),
    "lookup_gene": ToolSpec(_compute_lookup, ["genes", "enrichment", "data"]),
}


# --- execute ------------------------------------------------------------------


def execute(
    tool: str,
    inputs: dict[str, Path],
    parameters: dict[str, Any],
    upstream_runs: list[str] | None = None,
    use_cache: bool = True,
    runs_dir: Path = config.RUNS,
    results_root: Path = config.RESULTS,
) -> dict[str, Any]:
    """Run a registered tool with full provenance. Returns the run record.

    The record's "result" field holds the tool's structured output.
    """
    if tool not in TOOLS:
        raise KeyError(f"Unknown tool {tool!r}; registered tools: {sorted(TOOLS)}")
    spec = TOOLS[tool]
    t0 = time.perf_counter()
    parameters = to_jsonable(parameters)
    input_records = {k: {"path": rel(p), "sha256": data.sha256(Path(p))} for k, p in inputs.items()}
    software = software_versions()
    csha = code_sha(spec.modules)
    key = cache.cache_key(tool, {k: v["sha256"] for k, v in input_records.items()},
                          parameters, software, csha)

    stored = cache.load(tool, key, results_root) if use_cache else None
    cache_hit = stored is not None
    if not cache_hit:
        out_dir = cache.result_dir(tool, key, results_root)
        result = to_jsonable(spec.compute({k: Path(p) for k, p in inputs.items()}, parameters, out_dir))
        output_files = {k: rel(v) for k, v in result["artifacts"].items()
                        if isinstance(v, str)}
        stored = {
            "cache_key": key,
            "result": result,
            "output_files": output_files,
            "output_sha256": {k: data.sha256(absolute(v)) for k, v in output_files.items()},
        }
        cache.save(tool, key, stored, results_root)

    record = {
        "run_id": uuid.uuid4().hex[:8],
        "tool": tool,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "inputs": input_records,
        "parameters": parameters,
        "software": software,
        "code_sha256": csha,
        "git": git_state(),
        "outputs": stored["output_files"],
        "output_sha256": stored["output_sha256"],
        "summary": stored["result"].get("summary"),
        "upstream_runs": upstream_runs or [],
        "cache_key": key,
        "cache_hit": cache_hit,
        "duration_s": round(time.perf_counter() - t0, 3),
        "result": stored["result"],
    }
    runs_dir.mkdir(parents=True, exist_ok=True)
    (runs_dir / f"{record['run_id']}.json").write_text(json.dumps(record, indent=2))
    return record


def load_record(run_id: str, runs_dir: Path = config.RUNS) -> dict[str, Any]:
    """Load a run record, with a helpful error if it doesn't exist."""
    path = runs_dir / f"{run_id}.json"
    if not path.exists():
        known = sorted(p.stem for p in runs_dir.glob("*.json") if "." not in p.stem)
        raise FileNotFoundError(
            f"No run record {path}. Known run ids: {', '.join(known[-10:]) or 'none'}."
        )
    return json.loads(path.read_text())


# --- reproduce -----------------------------------------------------------------


def values_match(a: Any, b: Any, rtol: float = 1e-6, atol: float = 1e-9,
                 ignore: frozenset[str] = frozenset({"artifacts"})) -> list[str]:
    """Recursively compare two JSON-like values; return a list of differences."""
    if isinstance(a, dict) and isinstance(b, dict):
        diffs = []
        for k in sorted(set(a) | set(b)):
            if k in ignore:
                continue
            if k not in a or k not in b:
                diffs.append(f"key {k!r} only in {'original' if k in a else 'rerun'}")
            else:
                diffs += [f"{k}.{d}" for d in values_match(a[k], b[k], rtol, atol, ignore)]
        return diffs
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [f"length {len(a)} vs {len(b)}"]
        return [f"[{i}].{d}" for i, (x, y) in enumerate(zip(a, b))
                for d in values_match(x, y, rtol, atol, ignore)]
    if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        return [] if math.isclose(a, b, rel_tol=rtol, abs_tol=atol) else [f"{a} vs {b}"]
    return [] if a == b else [f"{a!r} vs {b!r}"]


def tables_match(a: Path, b: Path, rtol: float = 1e-6, atol: float = 1e-9) -> list[str]:
    """Compare two CSVs: identical columns/rows, numbers within tolerance."""
    x, y = pd.read_csv(a), pd.read_csv(b)
    if list(x.columns) != list(y.columns) or x.shape != y.shape:
        return [f"shape/columns differ: {x.shape} vs {y.shape}"]
    diffs = []
    for c in x.columns:
        if pd.api.types.is_numeric_dtype(x[c]) and pd.api.types.is_numeric_dtype(y[c]):
            ok = np.isclose(x[c], y[c], rtol=rtol, atol=atol, equal_nan=True)
        else:
            ok = (x[c].fillna("<NA>").astype(str) == y[c].fillna("<NA>").astype(str)).to_numpy()
        if not ok.all():
            diffs.append(f"column {c!r}: {int((~ok).sum())} rows differ")
    return diffs


def reproduce(run_id: str, runs_dir: Path = config.RUNS) -> dict[str, Any]:
    """Re-execute a run from its record (no cache) and verify outputs match.

    Fails early if any input file's checksum changed. Software version
    differences are reported as warnings. Tables and JSON are compared
    numerically (rtol 1e-6); plots by checksum, reported but not required.
    """
    rec = load_record(run_id, runs_dir)
    report: dict[str, Any] = {"run_id": run_id, "tool": rec["tool"], "checks": [], "warnings": []}

    for name, inp in rec["inputs"].items():
        path = absolute(inp["path"])
        if not path.exists():
            report["ok"] = False
            report["error"] = f"Input {name} ({inp['path']}) no longer exists; cannot reproduce."
            return report
        actual = data.sha256(path)
        if actual != inp["sha256"]:
            report["ok"] = False
            report["error"] = (f"Input {name} ({inp['path']}) changed since the run "
                               f"(sha256 {inp['sha256'][:12]}... -> {actual[:12]}...).")
            return report

    now = software_versions()
    for pkg, v in rec["software"].items():
        if now.get(pkg) != v:
            report["warnings"].append(f"{pkg}: recorded {v}, now {now.get(pkg)}")
    if code_sha(TOOLS[rec["tool"]].modules) != rec["code_sha256"]:
        report["warnings"].append("Tool source code changed since this run.")

    with tempfile.TemporaryDirectory() as tmp:
        t0 = time.perf_counter()
        spec = TOOLS[rec["tool"]]
        inputs = {k: absolute(v["path"]) for k, v in rec["inputs"].items()}
        rerun = to_jsonable(spec.compute(inputs, rec["parameters"], Path(tmp)))
        report["duration_s"] = round(time.perf_counter() - t0, 3)

        diffs = values_match(rec["result"], rerun)
        report["checks"].append({"output": "result", "method": "json (rtol 1e-6)",
                                 "match": not diffs, "detail": diffs[:5]})
        for name, path in rec["outputs"].items():
            new = Path(rerun["artifacts"][name])
            old = absolute(path)
            same_bytes = data.sha256(new) == rec["output_sha256"][name]
            if path.endswith(".csv"):
                d = [] if same_bytes else tables_match(old, new)
                report["checks"].append({"output": name, "method": "sha256" if same_bytes
                                         else "table (rtol 1e-6)", "match": not d, "detail": d[:5]})
            elif path.endswith(".json"):
                d = [] if same_bytes else values_match(json.loads(old.read_text()),
                                                       json.loads(new.read_text()))
                report["checks"].append({"output": name, "method": "sha256" if same_bytes
                                         else "json (rtol 1e-6)", "match": not d, "detail": d[:5]})
            else:
                report["checks"].append({"output": name, "method": "sha256 (informational)",
                                         "match": same_bytes, "required": False, "detail": []})

    report["ok"] = all(c["match"] for c in report["checks"] if c.get("required", True))
    report["timestamp_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    (runs_dir / f"{run_id}.reproduce.json").write_text(json.dumps(report, indent=2))
    return report


def log_event(kind: str, details: dict[str, Any], runs_dir: Path = config.RUNS) -> dict[str, Any]:
    """Record a non-computational event (e.g. a user overriding the QC gate) as a run record."""
    record = {
        "run_id": uuid.uuid4().hex[:8],
        "tool": kind,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "event": to_jsonable(details),
        "git": git_state(),
    }
    runs_dir.mkdir(parents=True, exist_ok=True)
    (runs_dir / f"{record['run_id']}.json").write_text(json.dumps(record, indent=2))
    return record
