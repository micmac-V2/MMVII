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
    # The sub-IFDs are assigned as plain dicts, not filled through get_ifd(): Exif.tobytes
    # walks only the top-level keys, so an IFD reached via get_ifd() alone is cached and
    # then dropped on save, leaving a JPEG with nothing but Make and Model.
    if sub:
        ex[0x8769] = dict(sub)
    if gps:
        ex[0x8825] = dict(gps)
    return ex


def _jpeg_xmp_segment(xmp: bytes) -> bytes:
    """An APP1 XMP segment, to be passed to Image.save as `extra`.

    Pillow's own `xmp=` save option is silently ignored for JPEG before 10.1, and
    pyproject allows Pillow 10.0, so the segment is built here instead: marker, a
    big-endian length that counts itself, the Adobe header, then the packet.
    """
    import struct

    payload = b"http://ns.adobe.com/xap/1.0/\x00" + xmp
    return b"\xff\xe1" + struct.pack(">H", len(payload) + 2) + payload


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
            root / f"DJI_20230421120000_{k + 1:04d}_D.JPG", exif=ex,
            extra=_jpeg_xmp_segment(xmp.encode()))
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


def make_canon_like_tiff(root: Path, name="IMG_9001.TIF", make="Canon", model="Canon EOS M",
                         focal=22.0, xres_per_inch=5798.66, size=(64, 48)):
    """A TIFF carrying the EXIF tags a Canon CR2 carries, for the raw IFD reader.

    CR2 is a TIFF container, so a plain TIFF with the same tags in the same places
    exercises `exif.read_tiff_tags` without needing a real raw file in the tree.
    """
    root.mkdir(parents=True, exist_ok=True)
    info = ImageFileDirectory_v2()
    info[0x010F] = make
    info[0x0110] = model
    info[0x8769] = {
        0x9003: "2026:11:10 15:01:16",
        0x920A: _r(focal),
        0xA20E: R(int(round(xres_per_inch * 100)), 100),  # FocalPlaneXResolution
        0xA210: 2,                                        # ... per inch
    }
    dest = root / name
    # Pillow picks its writer from the suffix and knows nothing of ".cr2", so write a
    # TIFF and then give it the requested name; the bytes are a TIFF container either way.
    tmp = dest.with_suffix(".tif")
    Image.new("I;16", size, 2048).save(tmp, tiffinfo=info)
    if tmp != dest:
        tmp.replace(dest)
    return dest


def bayer_mosaic(width=8, height=6, base=1000):
    """A small RGGB-ish mosaic whose four planes hold distinguishable constant-ish values.

    Plane (dy,dx) gets base*(1+index) plus the pixel's position, so a mis-strided split
    shows up immediately as a wrong value rather than a plausible one.
    """
    import numpy as np

    m = np.zeros((height, width), dtype=np.uint16)
    for dy in (0, 1):
        for dx in (0, 1):
            idx = dy * 2 + dx
            block = m[dy::2, dx::2]
            rows, cols = block.shape
            m[dy::2, dx::2] = (
                base * (1 + idx) + np.arange(rows * cols, dtype=np.uint16).reshape(rows, cols)
            )
    return m


# --- synthetic Blackmagic RAW container ------------------------------------------
#
# Enough of an ISO-BMFF/QuickTime file to exercise mmvii_studio.braw without a camera
# card: a video track for frame count and rate, two 20-byte IMU metadata tracks, and a
# moov/meta keys+ilst table of clip metadata.

def _box(typ: bytes, payload: bytes) -> bytes:
    import struct
    return struct.pack(">I", 8 + len(payload)) + typ + payload


def _imu_samples(values, tag: bytes) -> bytes:
    """Pack (N,3) floats as BRAW's 20-byte sample boxes."""
    import struct
    out = bytearray()
    for v in values:
        out += struct.pack(">I", 20) + tag + struct.pack("<3f", *v)
    return bytes(out)


def _stbl(fmt: bytes, sample_size: int, sample_count: int, timescale_delta: int,
          chunk_offset: int) -> bytes:
    import struct
    stsd = _box(b"stsd", struct.pack(">II", 0, 1) + _box(fmt, b""))
    stts = _box(b"stts", struct.pack(">III", 0, 1, sample_count) + struct.pack(">I", timescale_delta))
    stsz = _box(b"stsz", struct.pack(">III", 0, sample_size, sample_count))
    stsc = _box(b"stsc", struct.pack(">II", 0, 1) + struct.pack(">III", 1, sample_count, 1))
    stco = _box(b"stco", struct.pack(">II", 0, 1) + struct.pack(">I", chunk_offset))
    return _box(b"stbl", stsd + stts + stsz + stsc + stco)


def _trak(handler: bytes, fmt: bytes, timescale: int, duration: int, sample_size: int,
          sample_count: int, delta: int, chunk_offset: int) -> bytes:
    import struct
    mdhd = _box(b"mdhd", struct.pack(">IIIIIHH", 0, 0, 0, timescale, duration, 0, 0))
    hdlr = _box(b"hdlr", struct.pack(">I", 0) + b"mhlr" + handler + b"\0" * 12)
    minf = _box(b"minf", _stbl(fmt, sample_size, sample_count, delta, chunk_offset))
    return _box(b"trak", _box(b"mdia", mdhd + hdlr + minf))


def _keys_ilst(entries: list[tuple[str, int, bytes]]) -> bytes:
    """entries: [(key name, QuickTime value type, payload)]."""
    import struct
    keys = struct.pack(">II", 0, len(entries))
    for name, _vt, _pl in entries:
        nb = name.encode()
        keys += struct.pack(">I", 8 + len(nb)) + b"mdta" + nb
    ilst = b""
    for i, (_name, vt, pl) in enumerate(entries, start=1):
        data = _box(b"data", struct.pack(">II", vt, 0) + pl)
        ilst += struct.pack(">I", 8 + len(data)) + struct.pack(">I", i) + data
    hdlr = _box(b"hdlr", struct.pack(">I", 0) + b"\0" * 4 + b"mdta" + b"\0" * 12)
    return _box(b"meta", hdlr + _box(b"keys", keys) + _box(b"ilst", ilst))


def make_braw(path: Path, *, frames=48, fps_num=24, fps_den=1, imu_count=None,
              imu_rate=999, gyro=None, accel=None, camera_type="Blackmagic Pocket Cinema Camera 4K",
              firmware="8.1", width=1920, height=1080, line_time=7.246377,
              ois=1, lens="LUMIX G VARIO 14-42/F3.5-5.6", declared_imu_delta=1):
    """Write a synthetic .braw file. Returns the path.

    By default the IMU has one sample per 1/999 s spanning the clip, which is what a
    correct firmware writes. Pass `imu_count` smaller than that to reproduce the
    firmware 7.9.1 bug, where the track still declares 999 Hz.
    """
    import struct
    import numpy as np

    duration_s = frames * fps_den / fps_num
    if imu_count is None:
        imu_count = int(round(duration_s * imu_rate))
    if gyro is None:
        gyro = np.zeros((imu_count, 3), dtype=np.float32)
    if accel is None:
        # at rest, gravity along -Y as this camera's IMU reports it after negation
        accel = np.tile(np.array([0.0, -9.80665, 0.0], dtype=np.float32), (imu_count, 1))
    gyro = np.asarray(gyro, dtype=np.float32)[:imu_count]
    accel = np.asarray(accel, dtype=np.float32)[:imu_count]

    gyro_bytes = _imu_samples(gyro, b"mogy")
    # the file stores the pre-negation values, so flip back
    accel_bytes = _imu_samples(-accel, b"moac")
    video_bytes = b"\0" * (frames * 4)

    # mdat: video payload, then the two IMU chunks
    mdat_payload = video_bytes + gyro_bytes + accel_bytes
    # The file begins with the mdat box at offset 0, so its payload starts at byte 8.
    video_off = 8
    gyro_off = video_off + len(video_bytes)
    accel_off = gyro_off + len(gyro_bytes)
    mdat = _box(b"mdat", mdat_payload)

    traks = (
        _trak(b"vide", b"brst", fps_num, frames * fps_den, 4, frames, fps_den, video_off)
        + _trak(b"meta", b"mebx", 999, imu_count, 20, imu_count, declared_imu_delta, gyro_off)
        + _trak(b"meta", b"mebx", 999, imu_count, 20, imu_count, declared_imu_delta, accel_off)
    )
    meta = _keys_ilst([
        ("camera_type", 1, camera_type.encode()),
        ("firmware_version", 1, firmware.encode()),
        ("lens_type", 1, lens.encode()),
        ("crop_size", 71, struct.pack(">2f", width, height)),
        ("sensor_area_captured", 71, struct.pack(">2f", width, height)),
        ("crop_origin", 71, struct.pack(">2f", 16.0, 8.0)),
        ("sensor_line_time", 23, struct.pack(">f", line_time)),
        ("ois_enable", 76, struct.pack(">H", ois)),
        ("braw_compression_ratio", 1, b"8:1"),
        ("viewing_bmdgen", 76, struct.pack(">H", 5)),
    ])
    mvhd = _box(b"mvhd", struct.pack(">IIIII", 0, 0, 0, 1000, 1000) + b"\0" * 80)
    moov = _box(b"moov", mvhd + traks + meta)

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(mdat + moov)
    return path


def braw_frame_meta_blob(focal="14mm", aperture="f5.7", iso=1250, kelvin=5600,
                         shutter="172.8°", distance="135mm", gain=10.079369):
    """A bmdf per-frame metadata box, as it appears interleaved in mdat."""
    import struct

    def pad(s, n=24):
        b = s.encode("utf-8")
        return b + b"\0" * (n - len(b))

    inner = (
        _box(b"srte", struct.pack(">II", 24, 1))
        + _box(b"shtv", pad(shutter))
        + _box(b"agpf", struct.pack(">f", gain))
        + _box(b"aptr", pad(aperture))
        + _box(b"fcln", pad(focal))
        + _box(b"dsnc", pad(distance))
        + _box(b"isoe", struct.pack(">I", iso))
        + _box(b"wkel", struct.pack(">I", kelvin))
    )
    return _box(b"bmdf", inner)
