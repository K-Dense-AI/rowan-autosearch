from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


@pytest.fixture
def isolated_runs(tmp_path, monkeypatch):
    """Redirect all run state into a temp directory for each test."""
    from rowan_tools import state

    runs_dir = tmp_path / "runs"
    monkeypatch.setattr(state, "RUNS_DIR", runs_dir)
    return runs_dir
