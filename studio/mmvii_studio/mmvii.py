"""Locate and run the MMVII binary."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

_REPO_BIN = Path(__file__).resolve().parents[2] / "bin" / "MMVII"


def find_mmvii() -> str | None:
    """MMVII_BIN env var, then PATH, then the bin/ folder of this repository."""
    env = os.environ.get("MMVII_BIN")
    if env and Path(env).is_file():
        return env
    on_path = shutil.which("MMVII")
    if on_path:
        return on_path
    if _REPO_BIN.is_file():
        return str(_REPO_BIN)
    return None


class MMVIIError(RuntimeError):
    def __init__(self, cmd: list[str], output: str):
        super().__init__(f"MMVII failed: {' '.join(cmd)}")
        self.cmd = cmd
        self.output = output


def run(args: list[str], cwd: Path, timeout: float = 600) -> str:
    """Run `MMVII <args>` in cwd and return its combined output; raise MMVIIError on failure."""
    exe = find_mmvii()
    if exe is None:
        raise FileNotFoundError("MMVII binary not found (set MMVII_BIN or add it to PATH)")
    cmd = [exe, *args]
    proc = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True, timeout=timeout)
    if proc.returncode != 0:
        raise MMVIIError(cmd, proc.stdout)
    return proc.stdout
