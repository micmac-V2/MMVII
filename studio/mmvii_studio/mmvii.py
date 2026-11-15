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


# --- camera database and metadata rules -------------------------------------------
#
# MMVII needs a camera name, a focal length and a sensor geometry before it will build
# an internal calibration; without them OriCreateCalib stops with "Camera Name is not
# init". Images extracted from a Bayer mosaic carry no EXIF at all, so both have to be
# supplied explicitly: the geometry goes in a camera database, and a per-project rule
# attaches a model name and focal length to the files by filename pattern.
#
# Both are written to the LocalFolder database, which lives in the project's own
# MMVII-PhgrProj/MetaData/Std/ and is loaded last, so it overrides the global
# MMVII-RessourceDir/CameraDataBase.xml without modifying that shared file.

CAM_DB_LOCAL = "LocalFolder"

# eMTDIm values accepted by EditCalcMTDI. These are the *strings* MMVII parses, from the
# eMTDIm table in src/Serial/ElemStrToVal.cpp, not the C++ enum identifiers in
# include/MMVII_enums.h: all but one match once the leading "e" is dropped, but
# eNbPixel is spelled "NbPix" on the command line and in CalcMTD.xml.
MTD_MODES = ("Focalmm", "FocalEqui35", "FocalPix", "PPPix", "Aperture", "ModelCam",
             "NbPix", "AdditionalName")


def register_camera(
    project_dir: Path,
    name: str,
    pixel_size_micron: tuple[float, float],
    sensor_size_mm: tuple[float, float],
    image_size_px: tuple[int, int],
    *,
    db: str = CAM_DB_LOCAL,
    overwrite: bool = False,
) -> str:
    """Add a camera to an MMVII camera database. Returns MMVII's output.

    `overwrite` False adds only if the name is new (MMVII mode 1), True replaces an
    existing entry (mode 2).
    """
    params = [*pixel_size_micron, *sensor_size_mm, *image_size_px]
    param_arg = "Param=[" + ",".join(_fmt(v) for v in params) + "]"
    return run(
        ["EditCamDataBase", name, db, "2" if overwrite else "1", param_arg],
        cwd=project_dir,
    )


def set_metadata_rule(project_dir: Path, mode: str, value: str, pattern: str) -> str:
    """Attach `value` to every image matching `pattern`, for one metadata field.

    `mode` is one of :data:`MTD_MODES`; `pattern` is an MMVII regex matched against the
    image filename. The rule is appended to the project's CalcMTD.xml.
    """
    if mode not in MTD_MODES:
        raise ValueError(f"unknown metadata mode {mode!r}; expected one of {MTD_MODES}")
    return run(
        ["EditCalcMTDI", mode, f"Value={value}", f"PatternApply={pattern}", "Save=true"],
        cwd=project_dir,
    )


def register_plane_camera(project_dir: Path, cam, *, overwrite: bool = False) -> None:
    """Register a :class:`mmvii_studio.raw.PlaneCamera` and the rules that point at it.

    After this, ``MMVII OriCreateCalib "<pattern>" <Ori>`` succeeds on the plane images.
    """
    register_camera(
        project_dir,
        cam.name,
        cam.pixel_size_micron,
        cam.sensor_size_mm,
        cam.image_size_px,
        overwrite=overwrite,
    )
    set_metadata_rule(project_dir, "ModelCam", cam.name, cam.file_pattern)
    if cam.focal_mm is not None:
        set_metadata_rule(project_dir, "Focalmm", _fmt(cam.focal_mm), cam.file_pattern)


def _fmt(v) -> str:
    """Numbers for MMVII's command line: plain decimal, no exponent, no trailing noise."""
    if isinstance(v, int) or (isinstance(v, float) and v.is_integer()):
        return str(int(v))
    return f"{v:.6f}".rstrip("0").rstrip(".")
