"""Write small synthetic drone images with real EXIF/XMP, for tests and demos."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path

from PIL import Image
from PIL.TiffImagePlugin import IFDRational as R, ImageFileDirectory_v2

MICASENSE_XMP = ('<x:xmpmeta xmlns:x="adobe:ns:meta/"><rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">'
                 '<rdf:Description xmlns:Camera="http://pix4d.com/camera/1.0/" Camera:BandName="{band}" '
                 'Camera:RigCameraIndex="{idx}" Camera:CentralWavelength="{wl}"/></rdf:RDF></x:xmpmeta>')
DJI_XMP = ('<x:xmpmeta xmlns:x="adobe:ns:meta/"><rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">'
           '<rdf:Description xmlns:drone-dji="http://www.dji.com/drone-dji/1.0/" '
           'drone-dji:AbsoluteAltitude="+{abs_alt:.3f}" drone-dji:RelativeAltitude="+{rel_alt:.3f}"/>'
           '</rdf:RDF></x:xmpmeta>')
REDEDGE = [(1, "Blue", 475), (2, "Green", 560), (3, "Red", 668), (4, "NIR", 840), (5, "Red edge", 717)]


def _dms(v: float):
    v = abs(v)
    d = int(v)
    m = int((v - d) * 60)
    sec = Fraction((v - d - m / 60) * 3600).limit_denominator(100000)
    return (R(d), R(m), R(sec.numerator, sec.denominator))


def _r(v: float) -> R:
    f = Fraction(v).limit_denominator(1000)
    return R(f.numerator, f.denominator)


def _tags(make, model, lat, lon, alt, focal, f35=None, xres_per_cm=None):
    """(base tags, Exif sub-IFD, GPS sub-IFD) as plain dicts."""
    base = {0x010F: make, 0x0110: model}
    sub = {0x9003: "2023:04:21 12:00:00", 0x920A: _r(focal)}
    if f35:
        sub[0xA405] = f35
    if xres_per_cm:
        sub[0xA20E] = _r(xres_per_cm)
        sub[0xA210] = 3
    gps = {1: "N" if lat >= 0 else "S", 2: _dms(lat), 3: "E" if lon >= 0 else "W", 4: _dms(lon),
           5: b"\x00", 6: _r(alt)}
    return base, sub, gps


def _jpeg_exif(*args, **kw) -> Image.Exif:
    base, sub, gps = _tags(*args, **kw)
    ex = Image.Exif()
    ex.update(base)
    ex.get_ifd(0x8769).update(sub)
    ex.get_ifd(0x8825).update(gps)
    return ex


def _tiff_info(*args, xmp: bytes | None = None, **kw) -> ImageFileDirectory_v2:
    base, sub, gps = _tags(*args, **kw)
    info = ImageFileDirectory_v2()
    for k, v in base.items():
        info[k] = v
    info[0x8769] = sub
    info[0x8825] = gps
    if xmp:
        info[700] = xmp
    return info


def flight_positions(n: int, lat0=49.2500, lon0=-123.0000, step=0.0002):
    """A small serpentine grid of positions."""
    cols = max(2, int(n ** 0.5))
    out = []
    for k in range(n):
        r, c = divmod(k, cols)
        c = c if r % 2 == 0 else cols - 1 - c
        out.append((lat0 + r * step, lon0 + c * step))
    return out


def make_micasense(root: Path, n_captures=4, sub="0000SET/000"):
    folder = root / sub
    folder.mkdir(parents=True, exist_ok=True)
    for k, (lat, lon) in enumerate(flight_positions(n_captures)):
        for idx, band, wl in REDEDGE:
            img = Image.new("I;16", (64, 48), 1000 * idx)
            xmp = MICASENSE_XMP.format(band=band, idx=idx - 1, wl=wl).encode()
            info = _tiff_info("MicaSense", "RedEdge-M", lat, lon, 120.0 + k, 5.5, xres_per_cm=266.7, xmp=xmp)
            img.save(folder / f"IMG_{k:04d}_{idx}.tif", tiffinfo=info)
    return folder


def make_dji(root: Path, n=3):
    root.mkdir(parents=True, exist_ok=True)
    for k, (lat, lon) in enumerate(flight_positions(n, lat0=41.29, lon0=1.70)):
        ex = _jpeg_exif("DJI", "M3M", lat, lon, 250.0, 12.29, f35=24)
        xmp = DJI_XMP.format(abs_alt=262.5 + k, rel_alt=12.0)
        Image.new("RGB", (80, 60), (90, 120, 60)).save(
            root / f"DJI_20230421120000_{k + 1:04d}_D.JPG", exif=ex, xmp=xmp.encode())
        for tok in ("MS_G", "MS_NIR"):
            img = Image.new("I;16", (40, 30), 500)
            img.save(root / f"DJI_20230421120000_{k + 1:04d}_{tok}.TIF",
                     tiffinfo=_tiff_info("DJI", "M3M", lat, lon, 250.0, 4.34))
    return root


GCP_CSV = """Label,Easting,Northing,Elevation,Description
GCP 1,500010.250,5454600.120,12.345,corner
2,500040.000,5454630.500,13.000,
3,500070.125,5454610.875,12.800,fence
"""
