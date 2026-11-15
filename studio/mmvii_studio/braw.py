"""Read Blackmagic RAW (.braw) clip metadata and IMU motion data.

Why this exists
---------------
A Blackmagic Pocket Cinema Camera 4K (firmware 7.9 and later) records a real 6-DOF IMU
alongside the picture: a gyroscope in rad/s and an accelerometer in m/s². For
photogrammetry that is worth having twice over:

* the accelerometer sees gravity, which fixes roll and pitch in an absolute frame;
* ``sensor_line_time`` gives the rolling-shutter readout, which is a systematic
  geometric error in every frame of a moving CMOS camera.

Neither needs Blackmagic's SDK. A .braw file is an ISO-BMFF (QuickTime) container, so
the clip metadata sits in ``moov/meta`` as a ``keys``/``ilst`` pair and the IMU sits in
two timed-metadata tracks whose samples are 20-byte boxes: a big-endian size of 20, a
four-character type (``mogy`` for gyroscope, ``moac`` for accelerometer), then three
little-endian float32 values.

What this module does *not* do is decode pictures; that needs the Blackmagic SDK.

A firmware bug to know about
----------------------------
On firmware 7.9.1 the IMU tracks always declare a timescale of 999 with a sample delta
of 1, i.e. 999 Hz, whatever the real rate is. Measured against the video duration on a
Pocket 4K, the real rate was 999 Hz at 24.000 fps but 395.6 Hz at 23.976 fps, 664.1 Hz
at 29.97 fps and 1030.9 Hz at 60 fps. Firmware 8.1 clips were correct. Trusting the
declared rate on a 23.976 fps clip stretches IMU time by 2.5x, which ruins any
integration, so :func:`read_imu` derives the rate from the video duration by default and
reports both.
"""

from __future__ import annotations

import json
import struct
from dataclasses import dataclass, asdict, field
from pathlib import Path

import numpy as np

BRAW_SUFFIXES = {".braw"}

GYRO_TAG = b"mogy"
ACCEL_TAG = b"moac"
IMU_SAMPLE_SIZE = 20          # 4-byte size + 4-byte type + 3 float32
# One IMU sample is a complete little box: big-endian size, 4-char type, 3 LE float32.
IMU_DTYPE = np.dtype([("size", ">u4"), ("tag", "S4"), ("v", "<3f4")])
STANDARD_GRAVITY = 9.80665    # m/s², for sanity-checking the accelerometer

# QuickTime metadata value types seen in BRAW's ilst.
_QT_UTF8 = 1
_QT_FLOAT32 = 23
_QT_DIMENSIONS = 71           # two big-endian float32: width/height or x/y
_QT_INT = 76                  # Blackmagic writes small integers with this type
_QT_UINT32 = 77


def is_braw(path: str | Path) -> bool:
    """True when the suffix is .braw."""
    return Path(path).suffix.lower() in BRAW_SUFFIXES


# --- ISO-BMFF / QuickTime container walking --------------------------------------

def _boxes(f, start: int, end: int):
    """Yield (type, offset, size, header_size) for the boxes in [start, end)."""
    off = start
    while off + 8 <= end:
        f.seek(off)
        head = f.read(8)
        if len(head) < 8:
            return
        size, typ = struct.unpack(">I4s", head)
        hsize = 8
        if size == 1:                      # 64-bit extended size
            ext = f.read(8)
            if len(ext) < 8:
                return
            size = struct.unpack(">Q", ext)[0]
            hsize = 16
        elif size == 0:                    # runs to the end of its parent
            size = end - off
        if size < hsize or off + size > end:
            return
        yield typ.decode("latin1", "replace"), off, size, hsize
        off += size


def _find(f, start: int, end: int, want: str):
    """First child box of the given type, as (offset, size, header_size), or None."""
    for typ, off, size, hsize in _boxes(f, start, end):
        if typ == want:
            return off, size, hsize
    return None


def _full_box_entries(f, box, entry_struct: str):
    """Read a FullBox's (version/flags, count, entries) table."""
    off, size, hsize = box
    f.seek(off + hsize + 4)                # skip version/flags
    (count,) = struct.unpack(">I", f.read(4))
    width = struct.calcsize(entry_struct)
    limit = (off + size - (off + hsize + 8)) // max(width, 1)
    count = min(count, max(limit, 0))
    return [struct.unpack(entry_struct, f.read(width)) for _ in range(count)]


# --- clip metadata ---------------------------------------------------------------

@dataclass
class BrawClip:
    """The clip-level metadata of a .braw file."""

    path: str
    camera_type: str = ""
    camera_id: str = ""
    firmware_version: str = ""
    lens_type: str = ""
    manufacturer: str = ""
    clip_number: str = ""
    reel_name: str = ""
    date_recorded: str = ""
    compression: str = ""
    colour_science_gen: int | None = None
    width: int = 0
    height: int = 0
    crop_origin: tuple[float, float] | None = None
    frame_count: int = 0
    frame_rate: float | None = None
    sensor_line_time_us: float | None = None
    ois_enable: bool | None = None
    analog_gain: float | None = None

    @property
    def duration_s(self) -> float | None:
        if not self.frame_rate:
            return None
        return self.frame_count / self.frame_rate

    @property
    def readout_time_ms(self) -> float | None:
        """Rolling-shutter readout for the captured area: rows x line time.

        This is the time between the first and last sensor row of one frame. It is the
        quantity a rolling-shutter correction needs. Note that Gyroflow's parser halves
        it for a Pocket 4K on firmware 7.9, so treat the value as needing a check against
        a known motion before relying on it quantitatively.
        """
        if self.sensor_line_time_us is None or not self.height:
            return None
        return self.height * self.sensor_line_time_us / 1000.0

    def to_dict(self) -> dict:
        d = asdict(self)
        d["duration_s"] = self.duration_s
        d["readout_time_ms"] = self.readout_time_ms
        return d


def _read_keys_ilst(f, meta_box) -> dict:
    """Decode a QuickTime keys/ilst metadata pair into {key name: value}."""
    off, size, hsize = meta_box
    keys = _find(f, off + hsize, off + size, "keys")
    ilst = _find(f, off + hsize, off + size, "ilst")
    if not keys or not ilst:
        return {}

    ko, ksize, khsize = keys
    f.seek(ko + khsize + 4)                      # skip version/flags
    (count,) = struct.unpack(">I", f.read(4))
    names: list[str] = []
    for _ in range(min(count, 4096)):
        raw = f.read(4)
        if len(raw) < 4:
            break
        (klen,) = struct.unpack(">I", raw)
        if klen < 8:
            break
        f.read(4)                                # key namespace, always 'mdta' here
        names.append(f.read(klen - 8).decode("utf-8", "replace"))

    out: dict = {}
    io, isize, ihsize = ilst
    for typ, off2, size2, hsize2 in _boxes(f, io + ihsize, io + isize):
        # An ilst item's box type *is* its 1-based index into the keys table.
        index = struct.unpack(">I", typ.encode("latin1")[:4])[0]
        data = _find(f, off2 + hsize2, off2 + size2, "data")
        if not data:
            continue
        do, dsize, dhsize = data
        f.seek(do + dhsize)
        (vtype,) = struct.unpack(">I", f.read(4))
        f.read(4)                                # locale
        payload = f.read(dsize - dhsize - 8)
        name = names[index - 1] if 0 < index <= len(names) else f"#{index}"
        out[name] = _decode_qt_value(vtype, payload)
    return out


def _decode_qt_value(vtype: int, payload: bytes):
    if vtype == _QT_UTF8:
        return payload.decode("utf-8", "replace").strip()
    if vtype == _QT_FLOAT32 and len(payload) == 4:
        return struct.unpack(">f", payload)[0]
    if vtype == _QT_DIMENSIONS and len(payload) == 8:
        return struct.unpack(">2f", payload)
    if vtype in (_QT_INT, _QT_UINT32) and payload:
        return int.from_bytes(payload, "big")
    return payload


@dataclass
class _Track:
    handler: str
    format: str
    timescale: int
    duration: int
    sample_count: int
    sample_size: int
    stts: list[tuple[int, int]]
    chunk_offsets: list[int]
    stsc: list[tuple[int, int, int]]

    @property
    def declared_rate(self) -> float | None:
        if not self.stts or not self.stts[0][1]:
            return None
        return self.timescale / self.stts[0][1]


def _read_tracks(f, moov) -> list[_Track]:
    mo, msize, mhsize = moov
    tracks: list[_Track] = []
    for typ, off, size, hsize in _boxes(f, mo + mhsize, mo + msize):
        if typ != "trak":
            continue
        mdia = _find(f, off + hsize, off + size, "mdia")
        if not mdia:
            continue
        do, dsize, dhsize = mdia
        mdhd = _find(f, do + dhsize, do + dsize, "mdhd")
        hdlr = _find(f, do + dhsize, do + dsize, "hdlr")
        minf = _find(f, do + dhsize, do + dsize, "minf")
        if not (mdhd and hdlr and minf):
            continue
        f.seek(mdhd[0] + mdhd[2] + 12)
        timescale, duration = struct.unpack(">II", f.read(8))
        f.seek(hdlr[0] + hdlr[2] + 8)
        handler = f.read(4).decode("latin1", "replace")
        stbl = _find(f, minf[0] + minf[2], minf[0] + minf[1], "stbl")
        if not stbl:
            continue
        so, ssize, shsize = stbl

        fmt = ""
        stsd = _find(f, so + shsize, so + ssize, "stsd")
        if stsd:
            f.seek(stsd[0] + stsd[2] + 12)       # version/flags, entry count, entry size
            fmt = f.read(4).decode("latin1", "replace")

        sample_size = sample_count = 0
        stsz = _find(f, so + shsize, so + ssize, "stsz")
        if stsz:
            f.seek(stsz[0] + stsz[2] + 4)
            sample_size, sample_count = struct.unpack(">II", f.read(8))

        stts_box = _find(f, so + shsize, so + ssize, "stts")
        stts = _full_box_entries(f, stts_box, ">II") if stts_box else []
        stsc_box = _find(f, so + shsize, so + ssize, "stsc")
        stsc = _full_box_entries(f, stsc_box, ">III") if stsc_box else []

        offsets: list[int] = []
        stco = _find(f, so + shsize, so + ssize, "stco")
        co64 = _find(f, so + shsize, so + ssize, "co64")
        if co64:
            offsets = [e[0] for e in _full_box_entries(f, co64, ">Q")]
        elif stco:
            offsets = [e[0] for e in _full_box_entries(f, stco, ">I")]

        tracks.append(_Track(handler, fmt, timescale, duration, sample_count,
                             sample_size, stts, offsets, stsc))
    return tracks


def _chunk_plan(track: _Track) -> list[tuple[int, int]]:
    """[(chunk file offset, sample count)] for every chunk of a uniform-size track."""
    if not track.chunk_offsets or track.sample_size <= 0:
        return []
    entries = track.stsc or [(1, track.sample_count, 1)]
    per_chunk: list[int] = []
    for i, (first_chunk, samples, _desc) in enumerate(entries):
        last = entries[i + 1][0] - 1 if i + 1 < len(entries) else len(track.chunk_offsets)
        per_chunk.extend([samples] * max(0, last - first_chunk + 1))

    plan: list[tuple[int, int]] = []
    remaining = track.sample_count
    for ci, base in enumerate(track.chunk_offsets):
        if remaining <= 0:
            break
        n = per_chunk[ci] if ci < len(per_chunk) else (per_chunk[-1] if per_chunk else 0)
        n = min(n, remaining)
        if n <= 0:
            break
        plan.append((base, n))
        remaining -= n
    return plan


def _read_imu_track(f, track: _Track) -> tuple[bytes | None, np.ndarray]:
    """Read a uniform 20-byte metadata track as (tag, (N,3) float32).

    Samples inside one chunk are contiguous, so each chunk is a single read. That
    matters on a camera card, where thousands of small seeks are slow.
    """
    plan = _chunk_plan(track)
    if not plan:
        return None, np.empty((0, 3), dtype=np.float32)

    blocks: list[np.ndarray] = []
    tag: bytes | None = None
    for base, n in plan:
        f.seek(base)
        raw = f.read(n * IMU_SAMPLE_SIZE)
        got = len(raw) // IMU_SAMPLE_SIZE
        if got == 0:
            break
        rec = np.frombuffer(raw, dtype=IMU_DTYPE, count=got)
        if tag is None and got:
            tag = bytes(rec["tag"][0])
        blocks.append(np.asarray(rec["v"], dtype=np.float32))
        if got < n:
            break
    if not blocks:
        return tag, np.empty((0, 3), dtype=np.float32)
    return tag, np.concatenate(blocks, axis=0)


def _sample_offsets(track: _Track) -> list[int]:
    """Absolute file offset of every sample, from the chunk tables.

    Uses stsc (samples per chunk) with stco/co64 (chunk offsets) and a uniform stsz.
    Only uniform sample sizes are supported, which is what the IMU tracks use.
    """
    if not track.chunk_offsets or track.sample_size <= 0:
        return []
    per_chunk: list[int] = []
    entries = track.stsc or [(1, track.sample_count, 1)]
    for i, (first_chunk, samples, _desc) in enumerate(entries):
        last_chunk = entries[i + 1][0] - 1 if i + 1 < len(entries) else len(track.chunk_offsets)
        for _ in range(first_chunk, last_chunk + 1):
            per_chunk.append(samples)
    offsets: list[int] = []
    for ci, base in enumerate(track.chunk_offsets):
        n = per_chunk[ci] if ci < len(per_chunk) else (per_chunk[-1] if per_chunk else 0)
        for s in range(n):
            offsets.append(base + s * track.sample_size)
            if len(offsets) >= track.sample_count:
                return offsets
    return offsets


def read_clip(path: str | Path) -> BrawClip:
    """Clip-level metadata: camera, lens, geometry, frame rate, rolling-shutter timing."""
    path = Path(path)
    clip = BrawClip(path=str(path))
    with open(path, "rb") as f:
        f.seek(0, 2)
        size = f.tell()
        moov = _find(f, 0, size, "moov")
        if not moov:
            raise ValueError(f"{path.name}: no moov box; not an ISO-BMFF/BRAW file")

        meta = _find(f, moov[0] + moov[2], moov[0] + moov[1], "meta")
        md = _read_keys_ilst(f, meta) if meta else {}

        clip.camera_type = str(md.get("camera_type", "") or "")
        clip.camera_id = str(md.get("camera_id", "") or "")
        clip.firmware_version = str(md.get("firmware_version", "") or "")
        clip.lens_type = str(md.get("lens_type", "") or "").strip()
        clip.manufacturer = str(md.get("manufacturer", "") or "")
        clip.clip_number = str(md.get("clip_number", "") or "")
        clip.reel_name = str(md.get("reel_name", "") or "")
        clip.date_recorded = str(md.get("date_recorded", "") or "")
        clip.compression = str(md.get("braw_compression_ratio", "") or "")
        gen = md.get("viewing_bmdgen")
        clip.colour_science_gen = int(gen) if isinstance(gen, int) else None
        lt = md.get("sensor_line_time")
        clip.sensor_line_time_us = float(lt) if isinstance(lt, (int, float)) else None
        ois = md.get("ois_enable")
        clip.ois_enable = bool(ois) if isinstance(ois, int) else None
        gain = md.get("analog_gain")
        clip.analog_gain = float(gain) if isinstance(gain, (int, float)) else None

        size_tuple = md.get("crop_size") or md.get("sensor_area_captured")
        if isinstance(size_tuple, tuple) and len(size_tuple) == 2:
            clip.width, clip.height = int(size_tuple[0]), int(size_tuple[1])
        origin = md.get("crop_origin")
        if isinstance(origin, tuple) and len(origin) == 2:
            clip.crop_origin = (float(origin[0]), float(origin[1]))

        for tr in _read_tracks(f, moov):
            if tr.handler == "vide":
                clip.frame_count = tr.sample_count
                clip.frame_rate = tr.declared_rate
                break
    return clip


# --- IMU -------------------------------------------------------------------------

@dataclass
class ImuStream:
    """One IMU stream: N samples of three axes."""

    kind: str                    # "gyroscope" or "accelerometer"
    unit: str                    # "rad/s" or "m/s^2"
    samples: np.ndarray          # (N, 3) float32
    declared_rate_hz: float | None
    effective_rate_hz: float | None

    @property
    def count(self) -> int:
        return int(self.samples.shape[0])

    @property
    def rate_hz(self) -> float | None:
        return self.effective_rate_hz or self.declared_rate_hz

    def times(self) -> np.ndarray:
        """Sample times in seconds from the start of the clip."""
        rate = self.rate_hz
        if not rate:
            return np.zeros(self.count, dtype=np.float64)
        return np.arange(self.count, dtype=np.float64) / rate


@dataclass
class BrawImu:
    """Both IMU streams of a clip, plus the clip they came from."""

    clip: BrawClip
    gyro: ImuStream | None = None
    accel: ImuStream | None = None
    rate_source: str = "declared"      # "declared" or "derived from video duration"
    warnings: list[str] = field(default_factory=list)

    @property
    def has_imu(self) -> bool:
        return self.gyro is not None and self.accel is not None

    def gravity_check(self) -> float | None:
        """Median accelerometer magnitude, which should be close to 9.81 m/s² at rest.

        A value far from gravity means the clip was accelerating throughout, or the
        parse is wrong, or the sensor needs a scale calibration.
        """
        if self.accel is None or self.accel.count == 0:
            return None
        return float(np.median(np.linalg.norm(self.accel.samples, axis=1)))


def read_imu(path: str | Path, *, rate: str = "auto", negate_accel: bool = True) -> BrawImu:
    """Read the gyroscope and accelerometer tracks of a .braw clip.

    rate
        "auto" (default) derives the sample rate from the video duration whenever that
        disagrees with the track's declared rate by more than 2%, which works around the
        firmware 7.9.1 bug described in the module docstring. "declared" always trusts
        the container.
    negate_accel
        True (default) negates all three accelerometer axes, matching Gyroflow's
        telemetry-parser so that results are comparable with it.
    """
    path = Path(path)
    clip = read_clip(path)
    imu = BrawImu(clip=clip)

    with open(path, "rb") as f:
        f.seek(0, 2)
        size = f.tell()
        moov = _find(f, 0, size, "moov")
        if not moov:
            raise ValueError(f"{path.name}: no moov box")
        tracks = [t for t in _read_tracks(f, moov) if t.handler == "meta"]

        duration = clip.duration_s
        for tr in tracks:
            if tr.sample_size != IMU_SAMPLE_SIZE:
                continue
            tag, vals = _read_imu_track(f, tr)
            if tag not in (GYRO_TAG, ACCEL_TAG) or not len(vals):
                continue

            declared = tr.declared_rate
            effective = declared
            if rate == "auto" and duration and len(vals):
                measured = len(vals) / duration
                if declared and abs(measured - declared) / declared > 0.02:
                    effective = measured
                    imu.rate_source = "derived from video duration"
                    imu.warnings.append(
                        f"{path.name}: {tag.decode()} track declares "
                        f"{declared:.1f} Hz but spans {len(vals)} samples over "
                        f"{duration:.3f} s of video = {measured:.1f} Hz; using the "
                        f"measured rate (firmware {clip.firmware_version or '?'})")

            if tag == GYRO_TAG:
                imu.gyro = ImuStream("gyroscope", "rad/s", vals, declared, effective)
            else:
                if negate_accel:
                    vals = -vals
                imu.accel = ImuStream("accelerometer", "m/s^2", vals, declared, effective)

    return imu


# --- attitude from the IMU --------------------------------------------------------
#
# What follows turns the two IMU streams into one orientation per video frame. The
# gyroscope alone drifts; the accelerometer alone is useless while the camera
# accelerates. A complementary filter uses the gyroscope for the short term and lets the
# measured gravity direction pull roll and pitch back over the long term.
#
# Be clear about what this can and cannot give you:
#   * roll and pitch are anchored to gravity, so they do not drift;
#   * yaw (heading) is NOT observable from a gyroscope and accelerometer alone. It
#     starts at zero and drifts. Fixing heading needs a magnetometer, GNSS motion, or
#     the photogrammetry itself;
#   * the result is in the IMU's own axes. Using it on the *camera* needs a boresight
#     rotation between IMU and optical frame, which has to be calibrated. The default
#     here is identity, which is almost certainly wrong by some fixed rotation.


def _quat_mul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return np.array([
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
    ])


def _quat_conj(q: np.ndarray) -> np.ndarray:
    return np.array([q[0], -q[1], -q[2], -q[3]])


def _quat_rotate(q: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Rotate v by q (q maps body axes to world axes)."""
    qv = np.array([0.0, v[0], v[1], v[2]])
    return _quat_mul(_quat_mul(q, qv), _quat_conj(q))[1:]


def _quat_from_rotvec(r: np.ndarray) -> np.ndarray:
    theta = float(np.linalg.norm(r))
    if theta < 1e-12:
        return np.array([1.0, 0.0, 0.0, 0.0])
    axis = r / theta
    half = theta / 2.0
    return np.concatenate(([np.cos(half)], axis * np.sin(half)))


def _quat_between(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Shortest rotation taking unit vector a to unit vector b."""
    d = float(np.dot(a, b))
    if d > 1.0 - 1e-12:
        return np.array([1.0, 0.0, 0.0, 0.0])
    if d < -1.0 + 1e-12:                      # antiparallel: any perpendicular axis
        axis = np.cross(a, np.array([1.0, 0.0, 0.0]))
        if np.linalg.norm(axis) < 1e-9:
            axis = np.cross(a, np.array([0.0, 1.0, 0.0]))
        axis /= np.linalg.norm(axis)
        return np.concatenate(([0.0], axis))
    axis = np.cross(a, b)
    q = np.concatenate(([1.0 + d], axis))
    return q / np.linalg.norm(q)


def _quat_to_euler_deg(q: np.ndarray) -> np.ndarray:
    """Roll, pitch, yaw in degrees, as an intrinsic Z-Y-X (yaw-pitch-roll) sequence."""
    w, x, y, z = q
    sinr = 2.0 * (w * x + y * z)
    cosr = 1.0 - 2.0 * (x * x + y * y)
    roll = np.arctan2(sinr, cosr)
    sinp = np.clip(2.0 * (w * y - z * x), -1.0, 1.0)
    pitch = np.arcsin(sinp)
    siny = 2.0 * (w * z + x * y)
    cosy = 1.0 - 2.0 * (y * y + z * z)
    yaw = np.arctan2(siny, cosy)
    return np.degrees(np.array([roll, pitch, yaw]))


WORLD_DOWN = np.array([0.0, 0.0, -1.0])


@dataclass
class Attitude:
    """Orientation over time, from the IMU of one clip.

    :attr:`quat` is the authoritative representation: it is continuous, free of
    singularities, and consecutive samples are kept in the same hemisphere.
    :attr:`euler_deg` is derived from it for human reading and for formats that demand
    angles; it wraps at +/-180 degrees, which is a property of the representation and
    not a discontinuity in the underlying rotation.
    """

    times: np.ndarray            # (N,) seconds from clip start
    quat: np.ndarray             # (N, 4) w,x,y,z, body -> world
    euler_deg: np.ndarray        # (N, 3) roll, pitch, yaw -- for reading, not maths
    gravity_residual_deg: float  # median angle between measured and predicted gravity
    rate_hz: float | None
    yaw_is_absolute: bool = False   # never true from gyro+accel alone

    @property
    def count(self) -> int:
        return int(self.times.shape[0])


def integrate_attitude(imu: BrawImu, *, accel_gain: float = 2.0,
                       static_window: int = 200, reference: str = "start") -> Attitude:
    """Fuse the gyroscope and accelerometer into an orientation per IMU sample.

    accel_gain
        How hard the measured gravity pulls the estimate back, in units of 1/s. Too low
        and gyroscope drift wins; too high and linear acceleration is mistaken for tilt.
        Measured on Pocket 4K clips, the median gravity tracking error went 21.9 deg at
        0.05, 1.3 deg at 1.0, 0.65 deg at 2.0 and 0.45 deg at 5.0, so the default of 2.0
        is a deliberate middle. Check :attr:`Attitude.gravity_residual_deg` on your own
        footage rather than assuming.
    static_window
        How many leading samples to average for the initial orientation. The clip should
        begin reasonably still for this to be meaningful.
    reference
        "start" (default) reports orientation relative to the first frame, so the clip
        begins at zero roll, pitch and yaw. This is what you want for photogrammetry
        initialisation, and it avoids reporting the fixed rotation of the IMU inside the
        camera body as if it were camera tilt: on this Pocket 4K gravity lies along the
        IMU's Y axis, which otherwise shows up as a constant ~90 deg roll.
        "gravity" instead reports absolute roll and pitch in a world frame with Z up,
        which is only meaningful once the IMU-to-camera boresight is known.
    """
    if not imu.has_imu:
        raise ValueError("clip has no IMU data")
    gyro, accel = imu.gyro, imu.accel
    n = min(gyro.count, accel.count)
    if n == 0:
        raise ValueError("clip has empty IMU streams")
    rate = gyro.rate_hz
    if not rate:
        raise ValueError("no usable IMU sample rate")
    dt = 1.0 / rate

    g = np.asarray(accel.samples[:n], dtype=np.float64)
    w = np.asarray(gyro.samples[:n], dtype=np.float64)

    norms = np.linalg.norm(g, axis=1, keepdims=True)
    g_dir = np.divide(g, norms, out=np.zeros_like(g), where=norms > 1e-6)

    # Initial orientation: whatever rotation takes the measured gravity direction at the
    # start of the clip onto world down. Yaw is unobservable, so it starts at zero.
    lead = g_dir[: max(1, min(static_window, n))]
    lead_mean = lead.mean(axis=0)
    if np.linalg.norm(lead_mean) < 1e-6:
        q = np.array([1.0, 0.0, 0.0, 0.0])
    else:
        lead_mean /= np.linalg.norm(lead_mean)
        q = _quat_between(lead_mean, WORLD_DOWN)

    quats = np.empty((n, 4), dtype=np.float64)
    residuals = np.empty(n, dtype=np.float64)

    for i in range(n):
        # Predict with the gyroscope.
        q = _quat_mul(q, _quat_from_rotvec(w[i] * dt))
        q /= np.linalg.norm(q)

        # Correct with gravity, but only while the accelerometer plausibly sees it.
        meas = g_dir[i]
        if np.linalg.norm(meas) > 1e-6:
            pred = _quat_rotate(_quat_conj(q), WORLD_DOWN)   # world down, in body axes
            pred /= max(np.linalg.norm(pred), 1e-12)
            err = np.cross(pred, meas)                       # small-angle error, body axes
            residuals[i] = np.degrees(np.arcsin(np.clip(np.linalg.norm(err), -1.0, 1.0)))
            accel_is_trustworthy = abs(norms[i, 0] - STANDARD_GRAVITY) < 2.0
            if accel_is_trustworthy:
                # correction angle this step = gain (per second) x error x dt
                q = _quat_mul(q, _quat_from_rotvec(err * accel_gain * dt))
                q /= np.linalg.norm(q)
        else:
            residuals[i] = np.nan
        quats[i] = q

    # q and -q are the same rotation. Keep consecutive quaternions in the same
    # hemisphere so the track is continuous and safe to interpolate or difference.
    for i in range(1, n):
        if float(np.dot(quats[i], quats[i - 1])) < 0.0:
            quats[i] = -quats[i]

    if reference == "start":
        # Re-express every orientation relative to the first one.
        q0_inv = _quat_conj(quats[0] / np.linalg.norm(quats[0]))
        quats = np.array([_quat_mul(q0_inv, qq) for qq in quats])
    elif reference != "gravity":
        raise ValueError(f"reference must be 'start' or 'gravity', not {reference!r}")

    euler = np.array([_quat_to_euler_deg(qq) for qq in quats])
    return Attitude(
        times=np.arange(n, dtype=np.float64) * dt,
        quat=quats,
        euler_deg=euler,
        gravity_residual_deg=float(np.nanmedian(residuals)) if n else float("nan"),
        rate_hz=rate,
    )


def frame_attitude(imu: BrawImu, attitude: Attitude | None = None, *,
                   centre_of_exposure: bool = True) -> Attitude:
    """Resample an :class:`Attitude` to one value per video frame.

    Frame *i* is sampled at ``(i + 0.5) / fps`` when `centre_of_exposure` is set, which
    puts the sample at the middle of the frame rather than its first row. With a rolling
    shutter there is no single correct instant for a whole frame; see
    :attr:`BrawClip.readout_time_ms` for the spread.
    """
    attitude = attitude or integrate_attitude(imu)
    fps = imu.clip.frame_rate
    nframes = imu.clip.frame_count
    if not fps or not nframes:
        raise ValueError("clip has no frame rate or frame count")

    offset = 0.5 if centre_of_exposure else 0.0
    t = (np.arange(nframes, dtype=np.float64) + offset) / fps
    t = np.clip(t, attitude.times[0], attitude.times[-1])
    idx = np.searchsorted(attitude.times, t).clip(0, attitude.count - 1)

    return Attitude(
        times=t,
        quat=attitude.quat[idx],
        euler_deg=attitude.euler_deg[idx],
        gravity_residual_deg=attitude.gravity_residual_deg,
        rate_hz=fps,
    )


# --- export ----------------------------------------------------------------------

def write_attitude_csv(att: Attitude, path: str | Path, *, clip: BrawClip | None = None) -> Path:
    """Per-frame attitude as CSV: frame, time, roll, pitch, yaw, quaternion."""
    path = Path(path)
    lines = ["# frame,time_s,roll_deg,pitch_deg,yaw_deg,qw,qx,qy,qz"]
    if clip is not None:
        lines.insert(0, f"# {clip.camera_type} fw {clip.firmware_version} "
                        f"{clip.width}x{clip.height} @ {clip.frame_rate} fps")
        lines.insert(1, "# yaw is relative to the clip start and drifts; "
                        "roll and pitch are anchored to gravity")
    for i in range(att.count):
        r, p, y = att.euler_deg[i]
        qw, qx, qy, qz = att.quat[i]
        lines.append(f"{i},{att.times[i]:.6f},{r:.4f},{p:.4f},{y:.4f},"
                     f"{qw:.8f},{qx:.8f},{qy:.8f},{qz:.8f}")
    path.write_text("\n".join(lines) + "\n")
    return path


def write_mmvii_ori(att: Attitude, path: str | Path, *, name_template: str,
                    clip: BrawClip | None = None) -> Path:
    """Write an orientation file for ``MMVII ImportOri`` with format ``NXYZWPK``.

    Import it with, for example::

        MMVII ImportOri frames.ori NXYZWPK Calib InitOri AngU=degree

    **The positions are all zero.** An IMU measures orientation, not position: double
    integrating its accelerometer drifts hopelessly. MMVII's format demands X, Y and Z,
    so they are written as zeros and must be replaced from a real position source
    (GNSS, GCP-based resection) or left for the bundle adjustment to solve. Treat this
    file as initial *attitude* only.

    The angles are written as roll, pitch and yaw in degrees in the IMU frame, in the
    W, P, K columns. They are not omega/phi/kappa in MMVII's camera frame until a
    boresight rotation between the IMU and the optical axis has been applied, so expect
    a constant offset.
    """
    path = Path(path)
    lines = [
        "# MMVII ImportOri, format NXYZWPK, angles in degrees",
        "# WARNING X Y Z are placeholders (0): an IMU gives attitude, not position",
        "# WARNING W P K are IMU-frame roll/pitch/yaw, not boresighted to the camera",
        "# WARNING yaw is relative to the clip start and drifts (no magnetometer)",
    ]
    if clip is not None:
        lines.append(f"# {clip.camera_type} fw {clip.firmware_version} "
                     f"{Path(clip.path).name}")
    for i in range(att.count):
        r, p, y = att.euler_deg[i]
        lines.append(f"{name_template % i} 0 0 0 {r:.6f} {p:.6f} {y:.6f}")
    path.write_text("\n".join(lines) + "\n")
    return path


def imu_report(imu: BrawImu) -> dict:
    """A JSON-able summary of a clip's IMU, for logging and sanity checks."""
    rep: dict = {"clip": imu.clip.to_dict(), "warnings": list(imu.warnings),
                 "rate_source": imu.rate_source}
    grav = imu.gravity_check()
    rep["gravity_magnitude_m_s2"] = grav
    if grav:
        rep["gravity_scale_error_pct"] = 100.0 * (grav - STANDARD_GRAVITY) / STANDARD_GRAVITY
    for stream in (imu.gyro, imu.accel):
        if stream is None:
            continue
        rep[stream.kind] = {
            "unit": stream.unit,
            "count": stream.count,
            "declared_rate_hz": stream.declared_rate_hz,
            "effective_rate_hz": stream.effective_rate_hz,
            "abs_max_per_axis": [float(v) for v in np.abs(stream.samples).max(axis=0)]
            if stream.count else [],
        }
    return rep


# --- per-frame metadata ----------------------------------------------------------
#
# Each frame carries a 'bmdf' box of small sub-boxes. Decoded from a Pocket 4K clip:
#
#   srte  two u32: frame rate numerator and denominator (24/1)
#   shtv  string: shutter angle, e.g. "172.8°"
#   agpf  float32 big-endian: analog gain for this frame, in dB
#   aptr  string: aperture, e.g. "f5.7"
#   fcln  string: focal length, e.g. "14mm"
#   dsnc  string: focus distance, e.g. "135mm"
#   isoe  u32: ISO
#   wkel  u32: white balance in kelvin
#   wtin  u16: white balance tint
#
# Focal length is the one photogrammetry cannot do without, and on a zoom lens it can
# change during a clip, so it is worth reading per frame rather than once.

FRAME_META_BOX = b"bmdf"

_FRAME_STR_KEYS = {"shtv": "shutter_angle", "aptr": "aperture", "fcln": "focal_length",
                   "dsnc": "focus_distance"}
_FRAME_U32_KEYS = {"isoe": "iso", "wkel": "white_balance_kelvin"}

# Sensor dimensions are NOT in the file. These come from the manufacturer's published
# specifications and are needed to turn a focal length in mm into one in pixels. Verify
# against your own body before trusting a calibration built on them.
KNOWN_SENSORS: dict[str, tuple[float, float, int, int]] = {
    # camera_type: (sensor width mm, sensor height mm, native width px, native height px)
    "Blackmagic Pocket Cinema Camera 4K": (18.96, 10.00, 4096, 2160),
    "Blackmagic Pocket Cinema Camera 6K": (23.10, 12.99, 6144, 3456),
    "Blackmagic Pocket Cinema Camera 6K Pro": (23.10, 12.99, 6144, 3456),
    "Blackmagic Pocket Cinema Camera 6K G2": (23.10, 12.99, 6144, 3456),
}


def _parse_leading_number(text: str) -> float | None:
    """Pull the number off the front of strings like "14mm", "f5.7" or "172.8°"."""
    import re
    m = re.search(r"(\d+(?:\.\d+)?)", text)
    return float(m.group(1)) if m else None


def _decode_frame_meta(block: bytes) -> dict:
    """Decode one bmdf box payload into a plain dict."""
    out: dict = {}
    off = 0
    end = len(block)
    while off + 8 <= end:
        (bsize,) = struct.unpack_from(">I", block, off)
        btype = block[off + 4:off + 8].decode("latin1", "replace")
        if bsize < 8 or off + bsize > end:
            break
        payload = block[off + 8:off + bsize]
        if btype in _FRAME_STR_KEYS:
            text = payload.split(b"\x00")[0].decode("utf-8", "replace").strip()
            name = _FRAME_STR_KEYS[btype]
            out[name] = text
            num = _parse_leading_number(text)
            if num is not None:
                out[name + "_value"] = num
        elif btype in _FRAME_U32_KEYS and len(payload) >= 4:
            out[_FRAME_U32_KEYS[btype]] = struct.unpack_from(">I", payload, 0)[0]
        elif btype == "agpf" and len(payload) >= 4:
            out["analog_gain_db"] = round(struct.unpack_from(">f", payload, 0)[0], 6)
        elif btype == "srte" and len(payload) >= 8:
            num, den = struct.unpack_from(">II", payload, 0)
            out["frame_rate"] = (num / den) if den else None
        off += bsize
    return out


def read_frame_metadata(path: str | Path, *, max_frames: int | None = None,
                        max_bytes: int | None = None) -> list[dict]:
    """Per-frame metadata (focal length, aperture, ISO, shutter angle, white balance).

    The bmdf boxes are interleaved with the picture data in mdat, so this scans the file
    rather than using a sample table. Pass `max_frames` to stop early, which is what you
    want when you only need the focal length of a clip that does not zoom.
    """
    path = Path(path)
    frames: list[dict] = []
    chunk = 8 << 20
    overlap = 512                      # a bmdf box is ~256 bytes; never split one
    with open(path, "rb") as f:
        f.seek(0, 2)
        size = f.tell()
        limit = min(size, max_bytes) if max_bytes else size
        f.seek(0)
        prev = b""
        base = 0
        while base < limit:
            f.seek(base)
            buf = f.read(min(chunk, limit - base))
            if not buf:
                break
            hay = prev + buf
            start = 0
            while True:
                i = hay.find(FRAME_META_BOX, start)
                if i < 0 or i < 4:
                    break
                s = i - 4
                if s + 8 <= len(hay):
                    (bsize,) = struct.unpack_from(">I", hay, s)
                    if 8 < bsize <= 4096 and s + bsize <= len(hay):
                        meta = _decode_frame_meta(hay[s + 8:s + bsize])
                        if meta:
                            frames.append(meta)
                            if max_frames and len(frames) >= max_frames:
                                return frames
                start = i + 1
            prev = hay[-overlap:]
            base += len(buf)
    return frames


def pixel_size_mm(clip: BrawClip) -> float | None:
    """Sensor pixel pitch in mm, from the published sensor size for this camera.

    Returns None for a camera not in :data:`KNOWN_SENSORS`. The pitch is the *native*
    one: a windowed mode such as the Pocket 4K's 1920x1080 crops the sensor rather than
    downscaling it (its ``sensor_area_captured`` equals its ``crop_size``), so the pitch
    is unchanged by the mode.
    """
    spec = KNOWN_SENSORS.get(clip.camera_type)
    if not spec:
        return None
    width_mm, _h_mm, native_w, _native_h = spec
    if native_w <= 0:
        return None
    return width_mm / native_w


def focal_length_mm(path: str | Path) -> float | None:
    """Focal length of the first frame, in mm, or None when the lens did not report it."""
    frames = read_frame_metadata(path, max_frames=1, max_bytes=32 << 20)
    if not frames:
        return None
    return frames[0].get("focal_length_value")


def focal_length_px(clip: BrawClip, focal_mm: float | None) -> float | None:
    """Focal length in pixels, which is what an MMVII calibration holds."""
    pitch = pixel_size_mm(clip)
    if not pitch or not focal_mm:
        return None
    return focal_mm / pitch
