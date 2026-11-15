"""Tests for Blackmagic RAW ingest: container parsing, IMU, attitude, export."""

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))
import synth  # noqa: E402

from mmvii_studio import braw


def test_is_braw():
    assert braw.is_braw("A004_01010000_C041.braw")
    assert braw.is_braw(Path("x/y/CLIP.BRAW"))
    assert not braw.is_braw("clip.mov")
    assert not braw.is_braw("IMG_1753.CR2")


def test_read_clip(tmp_path):
    p = synth.make_braw(tmp_path / "c.braw", frames=48, width=1920, height=1080)
    c = braw.read_clip(p)
    assert c.camera_type == "Blackmagic Pocket Cinema Camera 4K"
    assert c.firmware_version == "8.1"
    assert c.lens_type.startswith("LUMIX")
    assert (c.width, c.height) == (1920, 1080)
    assert c.crop_origin == (16.0, 8.0)
    assert c.frame_count == 48
    assert c.frame_rate == pytest.approx(24.0)
    assert c.duration_s == pytest.approx(2.0)
    assert c.ois_enable is True
    assert c.colour_science_gen == 5
    assert c.compression == "8:1"
    assert c.sensor_line_time_us == pytest.approx(7.246377, rel=1e-5)


def test_readout_time_is_rows_times_line_time(tmp_path):
    p = synth.make_braw(tmp_path / "c.braw", height=1080, line_time=7.246377)
    c = braw.read_clip(p)
    assert c.readout_time_ms == pytest.approx(1080 * 7.246377 / 1000.0, rel=1e-6)


def test_read_clip_rejects_non_container(tmp_path):
    junk = tmp_path / "bad.braw"
    junk.write_bytes(b"not an ISO-BMFF file at all, no moov here")
    with pytest.raises(ValueError, match="moov"):
        braw.read_clip(junk)


def test_read_imu_units_and_negation(tmp_path):
    p = synth.make_braw(tmp_path / "c.braw", frames=24)
    imu = braw.read_imu(p)
    assert imu.has_imu
    assert imu.gyro.kind == "gyroscope" and imu.gyro.unit == "rad/s"
    assert imu.accel.kind == "accelerometer" and imu.accel.unit == "m/s^2"
    assert imu.gyro.count == imu.accel.count
    # synth writes gravity along -Y post-negation; the default read negates
    assert imu.accel.samples[0] == pytest.approx([0.0, -braw.STANDARD_GRAVITY, 0.0])
    assert imu.gravity_check() == pytest.approx(braw.STANDARD_GRAVITY, rel=1e-5)


def test_read_imu_without_negation(tmp_path):
    p = synth.make_braw(tmp_path / "c.braw", frames=24)
    imu = braw.read_imu(p, negate_accel=False)
    assert imu.accel.samples[0] == pytest.approx([0.0, braw.STANDARD_GRAVITY, 0.0])


def test_imu_times_use_the_rate(tmp_path):
    p = synth.make_braw(tmp_path / "c.braw", frames=24, imu_rate=999)
    imu = braw.read_imu(p)
    t = imu.gyro.times()
    assert t[0] == 0.0
    assert t[1] == pytest.approx(1 / 999.0)


def test_firmware_791_bug_rate_is_derived(tmp_path):
    """A short IMU track that still declares 999 Hz must not be trusted.

    This is the firmware 7.9.1 behaviour: on a 23.976 fps clip the track claimed 999 Hz
    while actually carrying about 396 Hz, which would stretch IMU time by 2.5x.
    """
    frames = 240
    p = synth.make_braw(tmp_path / "c.braw", frames=frames, fps_num=24000, fps_den=1001,
                        imu_count=3962, firmware="7.9.1")
    imu = braw.read_imu(p)
    duration = frames * 1001 / 24000
    assert imu.gyro.declared_rate_hz == pytest.approx(999.0)
    assert imu.gyro.effective_rate_hz == pytest.approx(3962 / duration, rel=1e-6)
    assert imu.gyro.effective_rate_hz < 500        # nowhere near the declared 999
    assert imu.rate_source == "derived from video duration"
    assert any("declares" in w for w in imu.warnings)


def test_declared_rate_can_be_forced(tmp_path):
    p = synth.make_braw(tmp_path / "c.braw", frames=240, fps_num=24000, fps_den=1001,
                        imu_count=3962, firmware="7.9.1")
    imu = braw.read_imu(p, rate="declared")
    assert imu.gyro.rate_hz == pytest.approx(999.0)
    assert imu.rate_source == "declared"
    assert imu.warnings == []


def test_correct_firmware_keeps_declared_rate(tmp_path):
    p = synth.make_braw(tmp_path / "c.braw", frames=48, imu_rate=999, firmware="8.1")
    imu = braw.read_imu(p)
    assert imu.rate_source == "declared"
    assert imu.warnings == []


def test_attitude_of_a_static_clip_is_zero(tmp_path):
    p = synth.make_braw(tmp_path / "c.braw", frames=48)
    att = braw.integrate_attitude(braw.read_imu(p))
    assert att.gravity_residual_deg == pytest.approx(0.0, abs=1e-6)
    assert np.allclose(att.euler_deg, 0.0, atol=1e-6)


def test_attitude_integrates_rotation_about_the_gravity_axis(tmp_path):
    """A known rotation must come back out.

    The rotation is about the gravity axis so the accelerometer correction has nothing
    to push against, which isolates the gyroscope integration.
    """
    rate, frames, fps = 999, 48, 24
    duration = frames / fps
    n = int(duration * rate)
    omega = 0.1                                   # rad/s about body Y, which is gravity
    gyro = np.tile(np.array([0.0, omega, 0.0], dtype=np.float32), (n, 1))
    p = synth.make_braw(tmp_path / "c.braw", frames=frames, imu_count=n, gyro=gyro)

    att = braw.integrate_attitude(braw.read_imu(p))
    total = 2 * math.degrees(math.acos(min(1.0, abs(float(att.quat[-1][0])))))
    expected = math.degrees(omega * (att.count - 1) / rate)
    assert total == pytest.approx(expected, rel=2e-3)
    # gravity is unchanged by this rotation, so the filter stays happy
    assert att.gravity_residual_deg == pytest.approx(0.0, abs=1e-3)


def test_attitude_quaternions_stay_in_one_hemisphere(tmp_path):
    n, frames = 4000, 96
    gyro = np.tile(np.array([0.0, 2.0, 0.0], dtype=np.float32), (n, 1))   # fast spin
    p = synth.make_braw(tmp_path / "c.braw", frames=frames, imu_count=n, gyro=gyro)
    att = braw.integrate_attitude(braw.read_imu(p))
    dots = np.sum(att.quat[1:] * att.quat[:-1], axis=1)
    assert (dots > 0).all()


def test_reference_start_vs_gravity(tmp_path):
    """"start" hides the IMU's mounting rotation; "gravity" exposes it."""
    p = synth.make_braw(tmp_path / "c.braw", frames=48)
    start = braw.integrate_attitude(braw.read_imu(p), reference="start")
    grav = braw.integrate_attitude(braw.read_imu(p), reference="gravity")
    assert np.allclose(start.euler_deg[0], 0.0, atol=1e-6)
    # synth puts gravity on -Y, so a Z-up world frame shows a ~90 deg roll
    assert abs(grav.euler_deg[0][0]) == pytest.approx(90.0, abs=1.0)


def test_reference_must_be_valid(tmp_path):
    p = synth.make_braw(tmp_path / "c.braw", frames=24)
    with pytest.raises(ValueError, match="reference"):
        braw.integrate_attitude(braw.read_imu(p), reference="nonsense")


def test_frame_attitude_has_one_sample_per_frame(tmp_path):
    p = synth.make_braw(tmp_path / "c.braw", frames=48)
    imu = braw.read_imu(p)
    fa = braw.frame_attitude(imu)
    assert fa.count == imu.clip.frame_count == 48
    assert fa.rate_hz == pytest.approx(24.0)
    # centre-of-exposure sampling puts frame 0 half a frame in
    assert fa.times[0] == pytest.approx(0.5 / 24.0, rel=1e-6)
    first_row = braw.frame_attitude(imu, centre_of_exposure=False)
    assert first_row.times[0] == pytest.approx(0.0)


def test_attitude_needs_imu(tmp_path):
    p = synth.make_braw(tmp_path / "c.braw", frames=24, imu_count=0)
    imu = braw.read_imu(p)
    with pytest.raises(ValueError):
        braw.integrate_attitude(imu)


def test_write_attitude_csv(tmp_path):
    p = synth.make_braw(tmp_path / "c.braw", frames=10)
    imu = braw.read_imu(p)
    fa = braw.frame_attitude(imu)
    out = braw.write_attitude_csv(fa, tmp_path / "att.csv", clip=imu.clip)
    lines = out.read_text().strip().splitlines()
    body = [l for l in lines if not l.startswith("#")]
    assert len(body) == 10
    assert lines[0].startswith("# Blackmagic")
    assert len(body[0].split(",")) == 9          # frame,time,r,p,y,qw,qx,qy,qz


def test_write_mmvii_ori_format_and_warnings(tmp_path):
    p = synth.make_braw(tmp_path / "c.braw", frames=5)
    imu = braw.read_imu(p)
    fa = braw.frame_attitude(imu)
    out = braw.write_mmvii_ori(fa, tmp_path / "f.ori", name_template="C041_%06d.tif",
                               clip=imu.clip)
    text = out.read_text()
    # the placeholder positions must be called out, not left to be discovered
    assert "placeholders" in text and "attitude, not position" in text
    assert "boresight" in text
    body = [l for l in text.splitlines() if not l.startswith("#")]
    assert len(body) == 5
    fields = body[0].split()
    assert len(fields) == 7                      # N X Y Z W P K
    assert fields[0] == "C041_000000.tif"
    assert fields[1:4] == ["0", "0", "0"]


def test_frame_metadata_decoding(tmp_path):
    blob = synth.braw_frame_meta_blob(focal="42mm", aperture="f8.0", iso=400,
                                      kelvin=5600, shutter="180.0°")
    f = tmp_path / "frames.braw"
    f.write_bytes(b"\0" * 64 + blob + b"\0" * 64 + blob)
    frames = braw.read_frame_metadata(f)
    assert len(frames) == 2
    m = frames[0]
    assert m["focal_length"] == "42mm" and m["focal_length_value"] == pytest.approx(42.0)
    assert m["aperture"] == "f8.0" and m["aperture_value"] == pytest.approx(8.0)
    assert m["iso"] == 400
    assert m["white_balance_kelvin"] == 5600
    assert m["shutter_angle_value"] == pytest.approx(180.0)
    assert m["frame_rate"] == pytest.approx(24.0)
    assert m["analog_gain_db"] == pytest.approx(10.079369, rel=1e-5)


def test_frame_metadata_max_frames(tmp_path):
    blob = synth.braw_frame_meta_blob()
    f = tmp_path / "frames.braw"
    f.write_bytes(blob * 5)
    assert len(braw.read_frame_metadata(f, max_frames=2)) == 2


def test_pixel_size_and_focal_in_pixels(tmp_path):
    p = synth.make_braw(tmp_path / "c.braw", frames=24, width=1920, height=1080)
    c = braw.read_clip(p)
    pitch = braw.pixel_size_mm(c)
    # Pocket 4K: 18.96 mm across 4096 native pixels, unchanged by a windowed mode
    assert pitch == pytest.approx(18.96 / 4096)
    assert braw.focal_length_px(c, 14.0) == pytest.approx(14.0 / pitch)
    assert braw.focal_length_px(c, None) is None


def test_pixel_size_unknown_camera(tmp_path):
    p = synth.make_braw(tmp_path / "c.braw", frames=24, camera_type="Some Other Camera")
    c = braw.read_clip(p)
    assert braw.pixel_size_mm(c) is None
    assert braw.focal_length_px(c, 14.0) is None


def test_imu_report_is_json_able(tmp_path):
    import json
    p = synth.make_braw(tmp_path / "c.braw", frames=24)
    rep = braw.imu_report(braw.read_imu(p))
    json.dumps(rep, default=str)                  # must not raise
    assert rep["gravity_magnitude_m_s2"] == pytest.approx(braw.STANDARD_GRAVITY, rel=1e-5)
    assert abs(rep["gravity_scale_error_pct"]) < 0.01
    assert rep["gyroscope"]["count"] == rep["accelerometer"]["count"]
