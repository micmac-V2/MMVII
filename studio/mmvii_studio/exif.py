"""Read the image metadata needed for project setup: GPS position, camera, focal, band.

Uses Pillow for EXIF and a small regex reader for XMP, which is where MicaSense and
DJI multispectral cameras store the band name and DJI stores its altitudes.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, asdict
from pathlib import Path

from PIL import Image

IMAGE_SUFFIXES = {".jpg", ".jpeg", ".tif", ".tiff", ".png"}

# EXIF tag ids
_TAG_MAKE = 0x010F
_TAG_MODEL = 0x0110
_IFD_EXIF = 0x8769
_IFD_GPS = 0x8825
_TAG_DATETIME_ORIG = 0x9003
_TAG_FOCAL = 0x920A
_TAG_FOCAL_35 = 0xA405
_TAG_FPLANE_XRES = 0xA20E
_TAG_FPLANE_UNIT = 0xA210
_TIFF_TAG_XMP = 700


@dataclass
class ImageInfo:
    path: str  # relative to the project root, with forward slashes
    width: int = 0
    height: int = 0
    make: str = ""
    model: str = ""
    datetime: str = ""
    focal_mm: float | None = None
    focal_35mm: float | None = None
    pixel_size_mm: float | None = None  # from focal plane resolution, if present
    lat: float | None = None
    lon: float | None = None
    alt: float | None = None  # GPS altitude (usually ellipsoidal or MSL, camera dependent)
    rel_alt: float | None = None  # DJI relative altitude above take-off
    band: str = ""  # band name from XMP, empty if none
    band_index: int | None = None  # 1-based band/rig index, if known

    @property
    def has_gps(self) -> bool:
        return self.lat is not None and self.lon is not None

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "ImageInfo":
        known = {k: v for k, v in d.items() if k in cls.__dataclass_fields__}
        return cls(**known)


def _ratio(v) -> float | None:
    if v is None:
        return None
    try:
        return float(v)
    except (TypeError, ValueError, ZeroDivisionError):
        try:
            num, den = v
            return float(num) / float(den) if den else None
        except Exception:
            return None


def _dms_to_deg(dms, ref) -> float | None:
    if not dms:
        return None
    parts = [_ratio(x) for x in dms]
    if any(p is None for p in parts):
        return None
    deg = parts[0] + (parts[1] if len(parts) > 1 else 0) / 60 + (parts[2] if len(parts) > 2 else 0) / 3600
    if isinstance(ref, bytes):
        ref = ref.decode(errors="ignore")
    if ref and ref.strip().upper() in ("S", "W"):
        deg = -deg
    return deg


def _clean(s) -> str:
    if isinstance(s, bytes):
        s = s.decode(errors="ignore")
    return str(s).strip("\x00 ").strip() if s is not None else ""


def parse_xmp(xmp: bytes | str | None) -> dict[str, str]:
    """Return XMP properties as {'ns:Key': value}, from both attribute and element forms."""
    if not xmp:
        return {}
    if isinstance(xmp, bytes):
        xmp = xmp.decode("utf-8", errors="ignore")
    out: dict[str, str] = {}
    for m in re.finditer(r'([A-Za-z][\w-]*:[A-Za-z][\w-]*)\s*=\s*"([^"]*)"', xmp):
        out.setdefault(m.group(1), m.group(2))
    for m in re.finditer(r"<([A-Za-z][\w-]*:[A-Za-z][\w-]*)>([^<]*)</\1>", xmp):
        out.setdefault(m.group(1), m.group(2).strip())
    return out


def _xmp_bytes(img: Image.Image):
    xmp = img.info.get("xmp") or img.info.get("XML:com.adobe.xmp")
    if xmp:
        return xmp
    tag = getattr(img, "tag_v2", None)
    if tag is not None:
        return tag.get(_TIFF_TAG_XMP)
    return None


def _float(s: str | None) -> float | None:
    if s is None or s == "":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def read_image_info(path: Path, root: Path) -> ImageInfo:
    rel = path.relative_to(root).as_posix()
    info = ImageInfo(path=rel)
    with Image.open(path) as img:
        info.width, info.height = img.size
        exif = img.getexif()
        info.make = _clean(exif.get(_TAG_MAKE))
        info.model = _clean(exif.get(_TAG_MODEL))

        sub = exif.get_ifd(_IFD_EXIF)
        info.datetime = _clean(sub.get(_TAG_DATETIME_ORIG))
        info.focal_mm = _ratio(sub.get(_TAG_FOCAL))
        f35 = _ratio(sub.get(_TAG_FOCAL_35))
        info.focal_35mm = f35 if f35 else None
        xres = _ratio(sub.get(_TAG_FPLANE_XRES))
        unit = sub.get(_TAG_FPLANE_UNIT)
        unit_mm = {2: 25.4, 3: 10.0, 4: 1.0}.get(unit)
        if xres and unit_mm:
            info.pixel_size_mm = unit_mm / xres

        gps = exif.get_ifd(_IFD_GPS)
        if gps:
            info.lat = _dms_to_deg(gps.get(2), gps.get(1))
            info.lon = _dms_to_deg(gps.get(4), gps.get(3))
            alt = _ratio(gps.get(6))
            if alt is not None:
                ref = gps.get(5)
                if ref in (1, b"\x01"):
                    alt = -alt
                info.alt = alt

        xmp = parse_xmp(_xmp_bytes(img))

    info.band = xmp.get("Camera:BandName", "")
    idx = xmp.get("Camera:RigCameraIndex")
    if idx is not None and idx.isdigit():
        info.band_index = int(idx) + 1
    abs_alt = _float(xmp.get("drone-dji:AbsoluteAltitude"))
    if abs_alt is not None:
        info.alt = abs_alt
    info.rel_alt = _float(xmp.get("drone-dji:RelativeAltitude"))
    return info


def find_images(root: Path) -> list[Path]:
    """All image files under root, skipping MMVII's own folders and hidden folders."""
    out = []
    for p in sorted(root.rglob("*")):
        if not p.is_file() or p.suffix.lower() not in IMAGE_SUFFIXES:
            continue
        parts = p.relative_to(root).parts
        if any(part.startswith(("MMVII-", ".")) for part in parts[:-1]):
            continue
        out.append(p)
    return out
