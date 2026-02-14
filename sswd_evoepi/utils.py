"""Utility functions for SSWD-EvoEpi.

General-purpose helpers: hashing, file I/O, timing, etc.
"""

from __future__ import annotations

import hashlib
import subprocess
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Generator


def file_sha256(path: str | Path) -> str:
    """Compute SHA-256 hex digest of a file."""
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            h.update(chunk)
    return h.hexdigest()


def config_hash(yaml_text: str) -> str:
    """SHA-256 of a YAML config string (for checkpoint tagging)."""
    return hashlib.sha256(yaml_text.encode('utf-8')).hexdigest()


def get_git_hash() -> str:
    """Return the current git commit hash, or 'unknown' if not in a repo."""
    try:
        result = subprocess.run(
            ['git', 'rev-parse', 'HEAD'],
            capture_output=True, text=True, timeout=5,
        )
        return result.stdout.strip() if result.returncode == 0 else 'unknown'
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return 'unknown'


@contextmanager
def timer(label: str = "") -> Generator[None, None, None]:
    """Simple context-manager timer. Prints elapsed time on exit."""
    start = time.perf_counter()
    yield
    elapsed = time.perf_counter() - start
    if label:
        print(f"[{label}] {elapsed:.3f}s")
    else:
        print(f"Elapsed: {elapsed:.3f}s")
