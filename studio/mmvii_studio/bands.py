"""Group the images of a multispectral camera into captures and bands.

A capture is one trigger of the camera; it holds one image per band. Supported naming:
  MicaSense:      IMG_0012_1.tif ... IMG_0012_5.tif     (band = rig index)
  DJI Mavic 3M:   DJI_20230609120000_0012_MS_NIR.TIF, ..._D.JPG (RGB)
Any other image becomes a single-band "RGB" capture, unless XMP gives a band name.
"""

from __future__ import annotations

import re
from collections import OrderedDict
from pathlib import PurePosixPath

from .exif import ImageInfo

# Band order of the MicaSense RedEdge / RedEdge-M (as in the Saczuk manual) and Altum.
MICASENSE_BANDS = {1: "Blue", 2: "Green", 3: "Red", 4: "NIR", 5: "RedEdge", 6: "LWIR"}
DJI_MS_BANDS = {"D": "RGB", "V": "RGB", "MS_G": "Green", "MS_R": "Red", "MS_RE": "RedEdge", "MS_NIR": "NIR"}

_RE_MICASENSE = re.compile(r"^(?P<cap>IMG_\d+)_(?P<band>\d+)$", re.I)
_RE_DJI = re.compile(r"^(?P<cap>DJI_\d+_\d+)_(?P<band>MS_[A-Z]+|D|V)$", re.I)

# Display order for band names
BAND_ORDER = ["RGB", "Blue", "Green", "Red", "RedEdge", "NIR", "LWIR", "Panchro"]


def band_sort_key(name: str):
    return (BAND_ORDER.index(name) if name in BAND_ORDER else len(BAND_ORDER), name)


def classify(info: ImageInfo) -> tuple[str, str]:
    """Return (capture_id, band_name) for one image."""
    p = PurePosixPath(info.path)
    folder = str(p.parent) if str(p.parent) != "." else ""
    stem = p.stem

    m = _RE_MICASENSE.match(stem)
    if m:
        idx = int(m.group("band"))
        band = info.band or MICASENSE_BANDS.get(idx, f"Band{idx}")
        if info.band_index is None:
            info.band_index = idx
        return _join(folder, m.group("cap")), _norm_band(band)

    m = _RE_DJI.match(stem)
    if m:
        token = m.group("band").upper()
        band = info.band or DJI_MS_BANDS.get(token, token)
        return _join(folder, m.group("cap")), _norm_band(band)

    return _join(folder, stem), _norm_band(info.band or "RGB")


def _join(folder: str, cap: str) -> str:
    return f"{folder}/{cap}" if folder else cap


def _norm_band(name: str) -> str:
    key = name.replace(" ", "").replace("-", "").replace("_", "").lower()
    aliases = {"rededge": "RedEdge", "re": "RedEdge", "nir": "NIR", "nearinfrared": "NIR",
               "blue": "Blue", "green": "Green", "red": "Red", "rgb": "RGB", "lwir": "LWIR",
               "panchro": "Panchro"}
    return aliases.get(key, name)


def group_captures(images: list[ImageInfo]) -> "OrderedDict[str, dict[str, ImageInfo]]":
    """capture_id -> {band: ImageInfo}, sorted by capture id."""
    caps: dict[str, dict[str, ImageInfo]] = {}
    for info in images:
        cap, band = classify(info)
        info.band = band
        caps.setdefault(cap, {})[band] = info
    return OrderedDict(sorted(caps.items()))


def bands_of(images: list[ImageInfo]) -> list[str]:
    return sorted({im.band or "RGB" for im in images}, key=band_sort_key)
