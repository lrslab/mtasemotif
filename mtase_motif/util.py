from __future__ import annotations

import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
import hashlib
from pathlib import Path
from typing import Optional, Sequence


class MtaseMotifError(RuntimeError):
    pass


@dataclass(frozen=True)
class CmdResult:
    argv: list[str]
    returncode: int
    stdout: str
    stderr: str


def which(tool: str) -> Optional[str]:
    env_tool = Path(sys.executable).resolve().parent / tool
    if env_tool.exists() and env_tool.is_file():
        return str(env_tool)
    return shutil.which(tool)


def run_cmd(
    argv: Sequence[str],
    *,
    cwd: Optional[Path] = None,
    check: bool = True,
    capture_output: bool = True,
    env: Optional[dict[str, str]] = None,
) -> CmdResult:
    if capture_output:
        proc = subprocess.run(
            list(argv),
            cwd=str(cwd) if cwd else None,
            capture_output=True,
            text=True,
            env=env or os.environ.copy(),
        )
        stdout = proc.stdout
        stderr = proc.stderr
    else:
        proc = subprocess.run(
            list(argv),
            cwd=str(cwd) if cwd else None,
            env=env or os.environ.copy(),
        )
        stdout = ""
        stderr = ""
    result = CmdResult(
        argv=list(argv),
        returncode=proc.returncode,
        stdout=stdout,
        stderr=stderr,
    )
    if check and proc.returncode != 0:
        raise MtaseMotifError(
            f"Command failed ({proc.returncode}): {' '.join(result.argv)}\n{result.stderr}"
        )
    return result


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()
