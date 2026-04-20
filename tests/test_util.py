from __future__ import annotations

import os
from pathlib import Path

from mtase_motif.util import which


def test_which_prefers_active_python_bin(monkeypatch, tmp_path: Path) -> None:
    env_bin = tmp_path / "env" / "bin"
    env_bin.mkdir(parents=True)
    env_tool = env_bin / "blastp"
    env_tool.write_text("#!/bin/sh\nexit 0\n")
    env_tool.chmod(0o755)

    monkeypatch.setattr("sys.executable", str(env_bin / "python"))
    monkeypatch.setenv("PATH", os.defpath)

    assert which("blastp") == str(env_tool)


def test_which_falls_back_to_path(monkeypatch, tmp_path: Path) -> None:
    env_bin = tmp_path / "env" / "bin"
    env_bin.mkdir(parents=True)

    path_bin = tmp_path / "pathbin"
    path_bin.mkdir()
    path_tool = path_bin / "hmmscan"
    path_tool.write_text("#!/bin/sh\nexit 0\n")
    path_tool.chmod(0o755)

    monkeypatch.setattr("sys.executable", str(env_bin / "python"))
    monkeypatch.setenv("PATH", str(path_bin))

    assert which("hmmscan") == str(path_tool)
