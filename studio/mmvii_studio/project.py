"""A Studio project: one folder holding the images, studio.json and MMVII-PhgrProj/.

MMVII is always run from the project folder, so image names are relative to it.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import warnings

from pyproj import CRS, Transformer
from pyproj.transformer import TransformerGroup
from pyproj.aoi import AreaOfInterest
from pyproj.database import query_utm_crs_info

from . import gcp as gcpmod
from .bands import bands_of, group_captures
from .exif import ImageInfo, find_images, read_image_info

PROJECT_FILE = "studio.json"
FORMAT_VERSION = 1


@dataclass
class Project:
    root: Path
    crs: str = ""  # project CRS, e.g. "EPSG:6653"; empty = automatic UTM
    images: list[ImageInfo] = field(default_factory=list)
    scan_errors: list[str] = field(default_factory=list)

    # ---------- persistence ----------
    @classmethod
    def open(cls, root: Path) -> "Project":
        root = Path(root).resolve()
        proj = cls(root=root)
        f = root / PROJECT_FILE
        if f.is_file():
            data = json.loads(f.read_text())
            proj.crs = data.get("crs", "")
            proj.images = [ImageInfo.from_dict(d) for d in data.get("images", [])]
            group_captures(proj.images)
        return proj

    def save(self):
        data = {"version": FORMAT_VERSION, "crs": self.crs,
                "images": [im.to_dict() for im in self.images]}
        (self.root / PROJECT_FILE).write_text(json.dumps(data, indent=1))

    # ---------- images ----------
    def scan_images(self, progress: Callable[[int, int, str], bool] | None = None):
        """Read metadata of every image under root. progress(i, n, name) returns False to cancel."""
        paths = find_images(self.root)
        images, errors = [], []
        for i, p in enumerate(paths):
            if progress and progress(i, len(paths), p.name) is False:
                break
            try:
                images.append(read_image_info(p, self.root))
            except Exception as e:  # unreadable or non-image file
                errors.append(f"{p.relative_to(self.root)}: {e}")
        group_captures(images)
        self.images, self.scan_errors = images, errors

    @property
    def bands(self) -> list[str]:
        return bands_of(self.images)

    @property
    def captures(self):
        return group_captures(self.images)

    def cameras_summary(self) -> list[dict]:
        """One row per distinct camera model / band / image size."""
        rows: dict[tuple, dict] = {}
        for im in self.images:
            key = (im.make, im.model, im.band, im.width, im.height, im.focal_mm)
            r = rows.setdefault(key, {"make": im.make, "model": im.model, "band": im.band,
                                      "width": im.width, "height": im.height, "focal_mm": im.focal_mm,
                                      "focal_35mm": im.focal_35mm, "pixel_size_mm": im.pixel_size_mm,
                                      "count": 0})
            r["count"] += 1
        return list(rows.values())

    # ---------- GCPs ----------
    def gcp_sets(self) -> list[gcpmod.GCPSet]:
        return gcpmod.list_mmvii_sets(self.root)

    # ---------- coordinate systems ----------
    def effective_crs(self) -> CRS:
        """The project CRS, or the UTM zone of the image GPS positions when unset."""
        if self.crs:
            return CRS.from_user_input(self.crs)
        gps = [im for im in self.images if im.has_gps]
        if not gps:
            return CRS.from_epsg(4326)
        lon = sum(im.lon for im in gps) / len(gps)
        lat = sum(im.lat for im in gps) / len(gps)
        return utm_for(lon, lat)


def utm_for(lon: float, lat: float) -> CRS:
    info = query_utm_crs_info(datum_name="WGS 84",
                              area_of_interest=AreaOfInterest(lon, lat, lon, lat))
    return CRS.from_epsg(int(info[0].code)) if info else CRS.from_epsg(4326)


def horizontal(crs: CRS) -> CRS:
    """The horizontal part of a compound CRS (e.g. NAD83(CSRS) UTM 10N of EPSG:6653)."""
    if crs.is_compound:
        for sub in crs.sub_crs_list:
            if not sub.is_vertical:
                return sub
    return crs


def to_map_xy(points_lonlat: list[tuple[float, float]], crs: CRS) -> list[tuple[float, float]]:
    tr = Transformer.from_crs(CRS.from_epsg(4326), horizontal(crs), always_xy=True)
    return [tr.transform(lon, lat) for lon, lat in points_lonlat]


def reproject_xy(points_xy: list[tuple[float, float]], src: CRS, dst: CRS) -> list[tuple[float, float]]:
    src_h, dst_h = horizontal(src), horizontal(dst)
    if src_h == dst_h:
        return list(points_xy)
    tr = Transformer.from_crs(src_h, dst_h, always_xy=True)
    return [tr.transform(x, y) for x, y in points_xy]


def missing_grids(src: CRS, dst: CRS) -> list[str]:
    """PROJ grids (e.g. a geoid model) needed for the best src->dst transformation but not installed.

    Without them PROJ falls back to a 'ballpark' transformation that, for heights, just copies
    ellipsoidal heights through. This is the Studio's version of the book's Geoids chapter.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tg = TransformerGroup(src, dst, always_xy=True)
    if tg.best_available or not tg.unavailable_operations:
        return []
    return [g.short_name for g in tg.unavailable_operations[0].grids if not g.available]
