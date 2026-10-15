"""Ground control points: read a delimited text file, import it with MMVII, read MMVII's result.

MMVII's ImportGCP splits fields on whitespace and commas, so a CSV with empty fields,
quoted values or spaces in names would shift columns. To avoid that we write a clean
"name x y z [sigma]" file into the project and import that one.
"""

from __future__ import annotations

import csv
import io
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

from pyproj import CRS
from pyproj.exceptions import CRSError

from . import mmvii

ROLES = ("name", "x", "y", "z", "sigma")
_GCP_DIR = Path("MMVII-PhgrProj") / "ObjCoordWorld"


def north_first(crs: str) -> bool:
    """True when the CRS's official axis order starts with northing/latitude (e.g. EPSG:4326).

    MMVII reads and writes coordinates in that official order, while the Studio always keeps
    x = easting/longitude and y = northing/latitude.
    """
    try:
        c = CRS.from_user_input(crs) if crs else None
    except CRSError:  # an MMVII alias such as L93 or RTL*...; PROJ cannot tell
        return False
    if c is None:
        return False
    if c.is_compound:
        c = next((s for s in c.sub_crs_list if not s.is_vertical), c)
    return bool(c.axis_info) and c.axis_info[0].direction.lower() in ("north", "south")


@dataclass
class GCP:
    name: str
    x: float  # easting / longitude
    y: float  # northing / latitude
    z: float
    sigma: float | None = None


@dataclass
class GCPSet:
    name: str
    crs: str
    points: list[GCP] = field(default_factory=list)


@dataclass
class TableSpec:
    """How to read a GCP text file. Column numbers are 0-based."""
    delimiter: str | None = None  # None = whitespace
    skip_rows: int = 0
    columns: dict[str, int] = field(default_factory=dict)  # role -> column


def sniff_delimiter(sample: str) -> str | None:
    lines = [ln for ln in sample.splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        return None
    for d in (",", ";", "\t"):
        if all(d in ln for ln in lines[: min(5, len(lines))]):
            return d
    return None


def read_rows(text: str, delimiter: str | None) -> list[list[str]]:
    rows = []
    if delimiter is None:
        for ln in text.splitlines():
            if ln.strip() and not ln.lstrip().startswith("#"):
                rows.append(ln.split())
    else:
        for row in csv.reader(io.StringIO(text), delimiter=delimiter):
            if any(c.strip() for c in row) and not row[0].lstrip().startswith("#"):
                rows.append([c.strip() for c in row])
    return rows


def guess_spec(text: str) -> TableSpec:
    """Guess delimiter, header rows and the name/x/y/z columns."""
    delim = sniff_delimiter(text)
    rows = read_rows(text, delim)
    skip = 0
    while skip < len(rows) and not _numeric_count(rows[skip]) >= 3:
        skip += 1
    spec = TableSpec(delimiter=delim, skip_rows=skip)
    if skip < len(rows):
        row = rows[skip]
        nums = [i for i, c in enumerate(row) if _is_num(c)]
        names = [i for i, c in enumerate(row) if not _is_num(c)]
        spec.columns["name"] = names[0] if names else 0
        coords = [i for i in nums if i != spec.columns["name"]]
        for role, col in zip(("x", "y", "z"), coords):
            spec.columns[role] = col
        if skip > 0:
            _apply_header_hints(spec, rows[skip - 1])
    return spec


_HEADER_HINTS = {
    "name": r"^(name|label|id|point|pt|gcp|marker)",
    "x": r"^(x|east|easting|e|lon|longitude)\b",
    "y": r"^(y|north|northing|n|lat|latitude)\b",
    "z": r"^(z|h|height|elev|elevation|alt|altitude)\b",
    "sigma": r"^(sigma|acc|accuracy|std)",
}


def _apply_header_hints(spec: TableSpec, header: list[str]):
    found = {}
    for i, h in enumerate(header):
        h = h.strip().lower()
        for role, pat in _HEADER_HINTS.items():
            if role not in found and re.match(pat, h):
                found[role] = i
    if {"name", "x", "y", "z"} <= found.keys():
        spec.columns = found


def _is_num(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def _numeric_count(row: list[str]) -> int:
    return sum(_is_num(c) for c in row)


def parse_points(text: str, spec: TableSpec) -> tuple[list[GCP], list[str]]:
    """Return points and a list of warnings (skipped lines, renamed points)."""
    missing = {"name", "x", "y", "z"} - spec.columns.keys()
    if missing:
        raise ValueError(f"Assign a column for: {', '.join(sorted(missing))}")
    pts, warnings = [], []
    for k, row in enumerate(read_rows(text, spec.delimiter)[spec.skip_rows:], start=spec.skip_rows + 1):
        try:
            name = row[spec.columns["name"]].strip()
            x, y, z = (float(row[spec.columns[r]]) for r in ("x", "y", "z"))
            sigma = float(row[spec.columns["sigma"]]) if "sigma" in spec.columns else None
        except (IndexError, ValueError):
            warnings.append(f"row {k}: could not read, skipped")
            continue
        clean = re.sub(r"[\s,]+", "_", name)
        if clean != name:
            warnings.append(f"row {k}: name '{name}' renamed to '{clean}'")
        pts.append(GCP(clean, x, y, z, sigma))
    return pts, warnings


def write_clean_file(points: list[GCP], path: Path, swap_xy: bool = False) -> str:
    """Write 'name x y z [sigma]' lines; return the MMVII format string for that file."""
    with_sigma = any(p.sigma is not None for p in points)
    with open(path, "w") as f:
        f.write("# written by MMVII Studio\n")
        for p in points:
            a, b = (p.y, p.x) if swap_xy else (p.x, p.y)
            line = f"{p.name} {a!r} {b!r} {p.z!r}"
            if with_sigma:
                line += f" {p.sigma if p.sigma is not None else 1.0!r}"
            f.write(line + "\n")
    return "NXYZSxyz" if with_sigma else "NXYZ"


def import_with_mmvii(root: Path, set_name: str, points: list[GCP], crs: str,
                      project_crs: str | None = None, sigma: float | None = None) -> str:
    """Create MMVII-PhgrProj/ObjCoordWorld/<set_name> from points. Returns MMVII's output."""
    if not re.fullmatch(r"[A-Za-z0-9_\-]+", set_name):
        raise ValueError("Set name may only contain letters, digits, '_' and '-'")
    src = root / f"GCP-{set_name}.txt"
    fmt = write_clean_file(points, src, swap_xy=north_first(crs))
    chsys = f"[{crs}]" if not project_crs or project_crs == crs else f"[{crs},{project_crs}]"
    args = ["ImportGCP", src.name, fmt, set_name, f"NameGCP={set_name}", f"ChSys={chsys}"]
    if sigma is not None and "Sxyz" not in fmt:
        args.append(f"Sigma={sigma}")
    return mmvii.run(args, cwd=root)


def _unquote(s: str | None) -> str:
    s = (s or "").strip()
    return s[1:-1] if len(s) >= 2 and s[0] == s[-1] == '"' else s


def read_mmvii_set(folder: Path) -> GCPSet:
    """Read an MMVII ObjCoordWorld folder (all MesGCP-*.xml files + CurSysCo.xml)."""
    crs = ""
    sysco = folder / "CurSysCo.xml"
    if sysco.is_file():
        crs = _unquote(ET.parse(sysco).getroot().findtext(".//SysCoData/Def"))
    gset = GCPSet(name=folder.name, crs=crs)
    for f in sorted(folder.glob("MesGCP-*.xml")):
        for el in ET.parse(f).getroot().iterfind(".//SetGCP/Measures/el"):
            name = _unquote(el.findtext("Name"))
            x, y, z = (float(v) for v in (el.findtext("Pt") or "").split()[:3])
            sig = (el.findtext("__Opt__Sigma2") or "").split()
            sigma = float(sig[0]) ** 0.5 if sig else None
            gset.points.append(GCP(name, x, y, z, sigma))
    if crs and north_first(crs):
        for g in gset.points:
            g.x, g.y = g.y, g.x
    return gset


def list_mmvii_sets(root: Path) -> list[GCPSet]:
    base = root / _GCP_DIR
    if not base.is_dir():
        return []
    return [read_mmvii_set(d) for d in sorted(base.iterdir()) if d.is_dir() and any(d.glob("MesGCP-*.xml"))]
