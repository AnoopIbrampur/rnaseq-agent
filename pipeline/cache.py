"""Content-addressed result cache.

A tool result is identified by a hash of everything that could change it:
the tool name, the SHA-256 of every input file, the parameters, the versions
of the scientific packages, and the SHA-256 of the tool's own source code.
Same key -> same result, so we can return the stored result instead of
recomputing. Results live in results/<tool>/<key[:12]>/ with a result.json.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from pipeline import config


def canonical_json(obj: Any) -> str:
    """Deterministic JSON (sorted keys, no whitespace) for hashing."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str)


def cache_key(
    tool: str,
    input_checksums: dict[str, str],
    parameters: dict[str, Any],
    software: dict[str, str],
    code_sha: str,
) -> str:
    """SHA-256 over the canonical JSON of everything that determines a result."""
    payload = {
        "tool": tool,
        "inputs": input_checksums,
        "parameters": parameters,
        "software": software,
        "code": code_sha,
    }
    return hashlib.sha256(canonical_json(payload).encode()).hexdigest()


def result_dir(tool: str, key: str, root: Path = config.RESULTS) -> Path:
    """Where outputs for this (tool, key) live."""
    return root / tool / key[:12]


def load(tool: str, key: str, root: Path = config.RESULTS) -> dict[str, Any] | None:
    """Return the stored result if present and all its output files still exist."""
    path = result_dir(tool, key, root) / "result.json"
    if not path.exists():
        return None
    stored = json.loads(path.read_text())
    if stored.get("cache_key") != key:
        return None
    if not all(Path(p).exists() for p in stored.get("output_files", {}).values()):
        return None
    return stored


def save(tool: str, key: str, stored: dict[str, Any], root: Path = config.RESULTS) -> Path:
    """Persist a result next to its output files."""
    path = result_dir(tool, key, root) / "result.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(stored, indent=2, default=str))
    return path
