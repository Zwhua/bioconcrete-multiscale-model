"""Machine-readable provenance manifests for formal analyses."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
from typing import Dict, Mapping, Optional, Sequence

from . import __version__
from .config import ModelConfig


def sha256_file(path: Path) -> str:
    """Return a streaming SHA-256 digest."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git(root: Path, *args: str) -> str:
    try:
        return subprocess.check_output(
            ["git", *args], cwd=str(root), text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def create_manifest(
    root: Path, config: ModelConfig, command: Sequence[str], random_seed: int,
    inputs: Optional[Mapping[str, Path]] = None, status: str = "running",
) -> Dict[str, object]:
    """Create a formal run manifest and warn when the worktree is dirty."""

    payload = json.dumps(config.to_dict(), sort_keys=True, separators=(",", ":"))
    dirty = bool(_git(root, "status", "--porcelain"))
    return {
        "model_version": __version__, "git_commit": _git(root, "rev-parse", "HEAD"),
        "git_worktree_dirty": dirty,
        "config_sha256": hashlib.sha256(payload.encode()).hexdigest(),
        "input_sha256": {name: sha256_file(path) for name, path in (inputs or {}).items()},
        "random_seed": int(random_seed), "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "completed_at_utc": None,
        "software_environment": {"python": sys.version.split()[0], "platform": platform.platform()},
        "command": list(command), "status": status,
        "release_warning": "dirty worktree; do not publish formal results" if dirty else "",
    }


def finish_manifest(manifest: Dict[str, object], status: str = "complete") -> Dict[str, object]:
    """Mark a manifest complete or failed."""

    result = dict(manifest)
    result["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
    result["status"] = status
    return result


def write_manifest(path: Path, manifest: Mapping[str, object]) -> None:
    """Write a stable JSON manifest."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(dict(manifest), indent=2), encoding="utf-8")
