"""Read camera raw files and split the Bayer mosaic into per-channel band images.

Why this exists
---------------
MMVII computes on single-channel images (``cIm2D``, ``eTN_U_INT2`` / ``eTN_REAL4``).
Its ``cRGBImage`` is three ``tU_INT1`` planes and is documented in
``include/MMVII_Image2D.h`` as "made for visu, not for intensive computation", so there
is no way to hand MMVII a colour image and have it use colour in matching.

Handing MMVII a raw file directly is worse than it looks: LibRaw-less readers (and MMVII
itself, through MMV1) fall back to the *embedded JPEG preview*, so a 14-bit CR2 arrives as
8-bit RGB and every extra level is discarded.

So instead of demosaicing to colour, we split the colour filter array into one
single-channel plane per filter. Each plane is half width and half height, carries the
sensor's full bit depth, and is written as a 16-bit TIFF that MMVII reads as ``U_INT2``.
The planes of one shot are named so that :mod:`mmvii_studio.bands` groups them into a
single capture with several bands, exactly as it does for a MicaSense rig.

No demosaicing happens anywhere here. Demosaicing interpolates between neighbouring
photosites, which invents detail that is correlated between channels and can bias tie
points; a split mosaic keeps every value a real measurement.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, asdict, field
from pathlib import Path

import numpy as np

# Raw suffixes we are willing to open. All are TIFF-ish container formats that LibRaw
# handles; the list is deliberately conservative rather than "everything LibRaw knows".
RAW_SUFFIXES = {
    ".cr2", ".cr3", ".crw",   # Canon
    ".nef", ".nrw",           # Nikon
    ".arw", ".srf", ".sr2",   # Sony
    ".raf",                   # Fujifilm (non-Bayer X-Trans: see is_bayer)
    ".orf",                   # Olympus
    ".rw2",                   # Panasonic
    ".pef",                   # Pentax
    ".srw",                   # Samsung
    ".dng",                   # Adobe / DJI / Leica
}

# Plane suffix -> band name used by bands.py. Two greens are kept apart on purpose:
# averaging them would be a one-line change but throws away a genuine second sample and
# the slight response difference between the two green photosites.
PLANE_BANDS = {"R": "Red", "G1": "Green1", "G2": "Green2", "B": "Blue"}


def is_raw(path: str | Path) -> bool:
    """True when the suffix is one of the raw formats handled here."""
    return Path(path).suffix.lower() in RAW_SUFFIXES


@dataclass
class RawInfo:
    """What the raw container says about the sensor, before any processing."""

    path: str
    raw_width: int = 0          # full sensor array, including masked borders
    raw_height: int = 0
    width: int = 0              # visible area
    height: int = 0
    bit_depth: int = 0          # inferred from white_level
    black_level: int = 0        # per-channel levels, when equal; else the mean
    black_level_per_channel: list[int] = field(default_factory=list)
    white_level: int = 0
    cfa_pattern: str = ""       # e.g. "RGGB", reading rows top-to-bottom
    planes: list[str] = field(default_factory=list)   # e.g. ["R","G1","G2","B"]
    camera_whitebalance: list[float] = field(default_factory=list)
    is_bayer: bool = True       # False for X-Trans and other non-2x2 mosaics

    @property
    def levels(self) -> int:
        """Distinct levels the container can hold above black."""
        return max(0, self.white_level - self.black_level)

    def to_dict(self) -> dict:
        return asdict(self)


def _plane_names(raw_pattern: np.ndarray, color_desc: bytes) -> tuple[str, list[str]]:
    """Map a 2x2 CFA index block to a pattern string and ordered plane names.

    ``raw_pattern`` holds colour indices into ``color_desc``; for a Canon RGGB sensor
    LibRaw reports ``[[0,1],[3,2]]`` with ``color_desc == b"RGBG"``, i.e. R G / G B.
    The two greens share the letter G, so they become G1 and G2 in raster order.
    """
    desc = color_desc.decode("ascii", "replace")
    letters = [desc[int(i)] if int(i) < len(desc) else "?" for i in raw_pattern.ravel()]
    pattern = "".join(letters)

    seen: dict[str, int] = {}
    names: list[str] = []
    counts = {ch: letters.count(ch) for ch in set(letters)}
    for ch in letters:
        if counts[ch] > 1:
            seen[ch] = seen.get(ch, 0) + 1
            names.append(f"{ch}{seen[ch]}")
        else:
            names.append(ch)
    return pattern, names


def read_raw_info(path: str | Path) -> RawInfo:
    """Open a raw file and report its sensor geometry and levels. Does not decode pixels."""
    import rawpy  # imported lazily so the GUI still starts without LibRaw

    path = Path(path)
    with rawpy.imread(str(path)) as r:
        sizes = r.sizes
        blacks = [int(b) for b in r.black_level_per_channel]
        white = int(r.white_level)
        pattern_arr = r.raw_pattern
        bayer = pattern_arr is not None and tuple(pattern_arr.shape) == (2, 2)
        if bayer:
            pattern, planes = _plane_names(pattern_arr, r.color_desc)
        else:
            pattern, planes = "", []
        return RawInfo(
            path=str(path),
            raw_width=int(sizes.raw_width),
            raw_height=int(sizes.raw_height),
            width=int(sizes.width),
            height=int(sizes.height),
            bit_depth=int(np.ceil(np.log2(white + 1))) if white > 0 else 0,
            black_level=int(round(sum(blacks) / len(blacks))) if blacks else 0,
            black_level_per_channel=blacks,
            white_level=white,
            cfa_pattern=pattern,
            planes=planes,
            camera_whitebalance=[float(x) for x in r.camera_whitebalance],
            is_bayer=bayer,
        )


def split_mosaic(
    mosaic: np.ndarray,
    raw_pattern: np.ndarray,
    color_desc: bytes,
    *,
    black_levels: list[int],
    white_level: int,
    normalize: bool = True,
    planes: list[str] | None = None,
) -> dict[str, np.ndarray]:
    """Split a 2-D CFA mosaic into one uint16 plane per colour filter.

    This is the core of the module with no LibRaw involved, so it can be exercised on a
    hand-built array. `raw_pattern` is the 2x2 block of colour indices into `color_desc`,
    as LibRaw reports them. Returns {plane name: array}, each half the mosaic's size.

    normalize
        True subtracts that plane's black level and scales the remainder to the full
        16-bit range, giving values linear in scene radiance. False returns the original
        sensor ADU.
    """
    if mosaic.ndim != 2:
        raise ValueError(f"expected a 2-D mosaic, got shape {mosaic.shape}")
    if tuple(raw_pattern.shape) != (2, 2):
        raise ValueError(f"expected a 2x2 CFA pattern, got {tuple(raw_pattern.shape)}")

    # Trim to an even number of rows/columns so every 2x2 block is whole.
    h, w = mosaic.shape
    mosaic = mosaic[: h - (h % 2), : w - (w % 2)]

    names = _plane_names(raw_pattern, color_desc)[1]
    wanted = set(names) if planes is None else set(planes)

    out: dict[str, np.ndarray] = {}
    for (dy, dx), name in zip([(0, 0), (0, 1), (1, 0), (1, 1)], names):
        if name not in wanted:
            continue
        plane = mosaic[dy::2, dx::2]
        # Per-plane black level, indexed by the CFA colour index at this position.
        cidx = int(raw_pattern[dy, dx])
        black = black_levels[cidx] if cidx < len(black_levels) else 0

        if normalize:
            span = white_level - black
            if span <= 0:
                raise ValueError(f"white level {white_level} <= black level {black}")
            # float32 represents every 14-bit integer exactly, so nothing is lost here.
            scaled = (plane.astype(np.float32) - black) * (65535.0 / span)
            out[name] = np.clip(scaled, 0, 65535).astype(np.uint16)
        else:
            out[name] = plane.astype(np.uint16, copy=True)
    return out


def _save_u16_tiff(arr: np.ndarray, path: Path, tags=None) -> None:
    """Write a 2-D uint16 array as a single-channel 16-bit TIFF (MMVII reads U_INT2)."""
    from PIL import Image

    if arr.dtype != np.uint16:
        raise ValueError(f"expected uint16, got {arr.dtype}")
    im = Image.fromarray(arr)
    if im.mode != "I;16":  # be explicit rather than trust the dtype mapping
        im = Image.frombytes("I;16", (arr.shape[1], arr.shape[0]), arr.tobytes())
    im.save(path, tiffinfo=tags) if tags is not None else im.save(path)


def plane_exif(path: Path, model_suffix: str = " Bayer"):
    """EXIF for an extracted plane, derived from the raw file's own EXIF.

    The plane is half the width and height of the sensor, so its pixels are twice as
    large and its focal-plane resolution is half the sensor's. Writing that into the
    file makes the plane self-describing: MMVII and Studio both read the camera model,
    the focal length and the pixel size straight from the TIFF, with no rule needed.

    The model gets a suffix so a plane is never mistaken for a full-frame image of the
    same camera, which has a different pixel count and pixel size.
    """
    from PIL.TiffImagePlugin import IFDRational, ImageFileDirectory_v2

    from .exif import read_tiff_tags  # local import: exif imports this module

    ifd0, exif, gps = read_tiff_tags(path)
    info = ImageFileDirectory_v2()

    make = ifd0.get(0x010F)
    model = ifd0.get(0x0110)
    if make:
        info[0x010F] = make
    if model:
        info[0x0110] = f"{model}{model_suffix}".strip()

    sub: dict = {}
    if 0x9003 in exif:
        sub[0x9003] = exif[0x9003]            # DateTimeOriginal
    if 0x920A in exif:
        num, den = exif[0x920A]
        sub[0x920A] = IFDRational(num, den)   # FocalLength, unchanged in mm
    if 0xA20E in exif and 0xA210 in exif:
        num, den = exif[0xA20E]
        # Half the resolution: one plane pixel spans two sensor pixels.
        sub[0xA20E] = IFDRational(num, den * 2)
        sub[0xA210] = exif[0xA210]
    if sub:
        info[0x8769] = sub

    # Only the position tags are worth carrying over. Copying the whole GPS IFD would
    # drag in GPSVersionID, which is a 4-byte array that Pillow refuses to write back.
    pos = {t: gps[t] for t in (1, 2, 3, 4, 5, 6) if t in gps}
    if 2 in pos and 4 in pos:
        info[0x8825] = pos
    return info


def extract_bayer_planes(
    path: str | Path,
    out_dir: str | Path,
    *,
    normalize: bool = True,
    planes: list[str] | None = None,
    overwrite: bool = False,
    write_exif: bool = True,
    model_suffix: str = " Bayer",
) -> list[Path]:
    """Split one raw file into per-channel 16-bit TIFFs and return the paths written.

    Output is ``<out_dir>/<stem>_<PLANE>.tif`` (``IMG_1753_R.tif`` and so on) plus a
    ``<stem>_raw.json`` sidecar recording the levels and geometry the planes came from.

    normalize
        True (default) subtracts the black level and scales the remaining range to the
        full 16-bit span, so the output is linear in scene radiance and uses the whole
        type. False keeps the original sensor ADU values, which is what you want when
        comparing against the camera's own numbers.
    planes
        Restrict to a subset, e.g. ``["G1","G2"]`` for the two green channels only.
    write_exif
        True (default) copies the camera, focal length and halved focal-plane resolution
        from the raw file into each plane, so the planes describe their own geometry.
    """
    import rawpy

    path = Path(path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    info = read_raw_info(path)
    if not info.is_bayer:
        raise ValueError(f"{path.name}: not a 2x2 Bayer mosaic (pattern {info.cfa_pattern!r})")

    wanted = set(planes) if planes else set(info.planes)
    unknown = wanted - set(info.planes)
    if unknown:
        raise ValueError(f"{path.name}: no such plane(s) {sorted(unknown)}; has {info.planes}")

    with rawpy.imread(str(path)) as r:
        planes_out = split_mosaic(
            r.raw_image_visible,
            np.asarray(r.raw_pattern),
            r.color_desc,
            black_levels=[int(b) for b in r.black_level_per_channel],
            white_level=int(r.white_level),
            normalize=normalize,
            planes=sorted(wanted),
        )

    tags = None
    if write_exif:
        try:
            tags = plane_exif(path, model_suffix)
        except Exception:
            tags = None  # no EXIF is survivable; the camera can be given by a rule

    written: list[Path] = []
    for name, plane in planes_out.items():
        dest = out_dir / f"{path.stem}_{name}.tif"
        if dest.exists() and not overwrite:
            written.append(dest)
            continue
        _save_u16_tiff(plane, dest, tags)
        written.append(dest)

    sidecar = out_dir / f"{path.stem}_raw.json"
    meta = info.to_dict()
    meta["normalized"] = normalize
    meta["plane_files"] = [p.name for p in written]
    meta["plane_bands"] = {n: PLANE_BANDS.get(n, n) for n in info.planes}
    sidecar.write_text(json.dumps(meta, indent=2) + "\n")

    return written


def plane_pixel_size_mm(sensor_pixel_size_mm: float | None) -> float | None:
    """Pixel size of an extracted plane: twice the sensor pitch, since the plane keeps
    one photosite per 2x2 block. MMVII needs this to turn a focal length in mm into
    a focal length in pixels."""
    if sensor_pixel_size_mm is None:
        return None
    return sensor_pixel_size_mm * 2.0


def plane_file_pattern(planes: list[str], suffix: str = "tif") -> str:
    """MMVII-style regex matching the plane files of `planes`.

    For the usual Canon planes this gives ``.*_(G1|G2|B|R)\\.tif``, which is what goes
    into an MMVII metadata rule's PatternApply.
    """
    alts = "|".join(sorted(planes, key=lambda n: (-len(n), n)))
    return r".*_(" + alts + r")\." + suffix.lstrip(".")


@dataclass
class PlaneCamera:
    """The camera geometry MMVII needs for images extracted from a Bayer mosaic.

    A plane holds one photosite per 2x2 block, so it has half the pixel count in each
    axis and twice the pixel pitch, over the same physical sensor area. Focal length in
    mm is unchanged; focal length in *pixels* therefore halves, which is why the plane
    needs its own camera entry rather than reusing the full-frame one.
    """

    name: str
    pixel_size_micron: tuple[float, float]
    sensor_size_mm: tuple[float, float]
    image_size_px: tuple[int, int]
    focal_mm: float | None
    file_pattern: str

    @property
    def focal_px(self) -> float | None:
        """Focal length in plane pixels, the quantity MMVII puts in a calibration."""
        if self.focal_mm is None or self.pixel_size_micron[0] <= 0:
            return None
        return self.focal_mm / (self.pixel_size_micron[0] / 1000.0)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["focal_px"] = self.focal_px
        return d


def bayer_camera(
    info: RawInfo,
    sensor_pixel_size_mm: float,
    model: str,
    focal_mm: float | None = None,
    suffix: str = " Bayer",
) -> PlaneCamera:
    """Camera geometry for the planes extracted from `info`.

    `sensor_pixel_size_mm` is the full-sensor pitch, normally from the EXIF focal-plane
    resolution; `model` is the EXIF camera model. The returned name is `model + suffix`
    so the plane camera cannot be confused with the full-frame one.
    """
    if sensor_pixel_size_mm <= 0:
        raise ValueError("sensor_pixel_size_mm must be positive")

    plane_pitch_mm = sensor_pixel_size_mm * 2.0
    # Plane size uses the same even-trim as extract_bayer_planes.
    pw = (info.width - info.width % 2) // 2
    ph = (info.height - info.height % 2) // 2
    return PlaneCamera(
        name=f"{model}{suffix}".strip(),
        pixel_size_micron=(plane_pitch_mm * 1000.0, plane_pitch_mm * 1000.0),
        sensor_size_mm=(pw * plane_pitch_mm, ph * plane_pitch_mm),
        image_size_px=(pw, ph),
        focal_mm=focal_mm,
        file_pattern=plane_file_pattern(info.planes or list(PLANE_BANDS)),
    )
