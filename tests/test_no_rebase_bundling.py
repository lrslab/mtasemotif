from __future__ import annotations

import subprocess
import sys
import tarfile
from pathlib import Path


def test_sdist_does_not_bundle_rebase_data(tmp_path: Path) -> None:
    subprocess.check_call([sys.executable, "setup.py", "sdist", "--dist-dir", str(tmp_path)])
    tgz = next(tmp_path.glob("*.tar.gz"))
    with tarfile.open(tgz, "r:gz") as tf:
        names = tf.getnames()

    forbidden = [n for n in names if "emboss_" in n or "/rebase/" in n]
    assert not forbidden, f"Forbidden REBASE-like files bundled in sdist: {forbidden[:20]}"
