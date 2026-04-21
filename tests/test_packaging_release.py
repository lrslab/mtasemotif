from __future__ import annotations

import os
import re
import subprocess
import sys
import tarfile
from pathlib import Path


def _build_sdist(tmp_path: Path) -> Path:
    subprocess.check_call([sys.executable, "-m", "build", "--sdist", "--outdir", str(tmp_path)])
    return next(tmp_path.glob("*.tar.gz"))


def test_sdist_excludes_cache_and_macos_archive_dirs(tmp_path: Path) -> None:
    tgz = _build_sdist(tmp_path)
    with tarfile.open(tgz, "r:gz") as tf:
        names = tf.getnames()

    forbidden = [
        name
        for name in names
        if "__pycache__" in name or "__MACOSX" in name or name.endswith(".pyc")
    ]
    assert not forbidden, f"Unexpected transient files bundled in sdist: {forbidden[:20]}"


def test_sdist_excludes_local_db_and_run_artifacts(tmp_path: Path) -> None:
    tgz = _build_sdist(tmp_path)
    with tarfile.open(tgz, "r:gz") as tf:
        names = tf.getnames()

    forbidden_fragments = (
        "/real_db/",
        "/real_db_sources/",
        "/real_run_realdb",
        "/benchmarks/known_motif_panel/runs/",
        "/_db_rebase_test_",
        "/_db_test_",
        "/_rebase_src_",
    )
    forbidden = [name for name in names if any(fragment in name for fragment in forbidden_fragments)]
    assert not forbidden, f"Unexpected local data bundled in sdist: {forbidden[:20]}"


def test_module_invocation_shows_help() -> None:
    env = dict(os.environ, NO_COLOR="1", TERM="dumb")
    proc = subprocess.run(
        [sys.executable, "-m", "mtase_motif.cli", "--help"],
        capture_output=True,
        text=True,
        env=env,
    )
    stdout = re.sub(r"\x1b\[[0-9;?]*[ -/]*[@-~]", "", proc.stdout)

    assert proc.returncode == 0
    assert "Usage:" in stdout
    assert "mtase_motif.cli" in stdout
    assert "run" in stdout
    assert "db" in stdout
