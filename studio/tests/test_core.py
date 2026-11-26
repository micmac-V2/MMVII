import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
import synth  # noqa: E402

from PIL import Image

from mmvii_studio import gcp as gcpmod
from mmvii_studio import mmvii
from mmvii_studio.bands import classify, group_captures
from mmvii_studio.exif import ImageInfo, parse_xmp, read_image_info
from mmvii_studio.project import Project, horizontal, to_map_xy
from pyproj import CRS


def test_micasense_exif_and_bands(tmp_path):
    synth.make_micasense(tmp_path, n_captures=3)
    proj = Project.open(tmp_path)
    proj.scan_images()
    assert len(proj.images) == 15
    assert proj.bands == ["Blue", "Green", "Red", "RedEdge", "NIR"]
    caps = proj.captures
    assert list(caps) == [f"0000SET/000/IMG_{k:04d}" for k in range(3)]
    im = caps["0000SET/000/IMG_0001"]["NIR"]
    assert im.band_index == 4
    assert im.make == "MicaSense" and im.model == "RedEdge-M"
    assert im.focal_mm == pytest.approx(5.5)
    assert im.pixel_size_mm == pytest.approx(10 / 266.7, rel=1e-3)
    lat, lon = synth.flight_positions(3)[1]
    assert im.lat == pytest.approx(lat, abs=1e-7) and im.lon == pytest.approx(lon, abs=1e-7)
    assert im.alt == pytest.approx(121.0)


def test_dji_xmp_altitude_and_bands(tmp_path):
    synth.make_dji(tmp_path, n=2)
    proj = Project.open(tmp_path)
    proj.scan_images()
    assert proj.bands == ["RGB", "Green", "NIR"]
    rgb = proj.captures["DJI_20230421120000_0002"]["RGB"]
    assert rgb.alt == pytest.approx(263.5)  # XMP AbsoluteAltitude wins over EXIF
    assert rgb.rel_alt == pytest.approx(12.0)
    assert rgb.focal_35mm == 24
    assert rgb.lon == pytest.approx(synth.flight_positions(2, 41.29, 1.70)[1][1], abs=1e-7)


def test_project_roundtrip_and_auto_utm(tmp_path):
    synth.make_micasense(tmp_path, n_captures=2)
    p = Project.open(tmp_path)
    p.scan_images()
    p.save()
    q = Project.open(tmp_path)
    assert [im.to_dict() for im in q.images] == [im.to_dict() for im in p.images]
    assert q.effective_crs().to_epsg() == 32610  # Vancouver area -> UTM 10N
    q.crs = "EPSG:6653"
    assert horizontal(q.effective_crs()).to_epsg() == 3157


def test_classify_fallback():
    assert classify(ImageInfo(path="flight/P0001.JPG")) == ("flight/P0001", "RGB")


def test_xmp_element_form():
    assert parse_xmp(b"<Camera:BandName>NIR</Camera:BandName>") == {"Camera:BandName": "NIR"}


def test_map_xy_utm():
    (x, y), = to_map_xy([(-123.0, 49.25)], CRS.from_epsg(32610))
    assert x == pytest.approx(500000, abs=1)  # central meridian of zone 10
    assert y == pytest.approx(5455000, abs=2000)


def test_gcp_guess_and_parse():
    spec = gcpmod.guess_spec(synth.GCP_CSV)
    assert spec.delimiter == "," and spec.skip_rows == 1
    assert spec.columns == {"name": 0, "x": 1, "y": 2, "z": 3}
    pts, warns = gcpmod.parse_points(synth.GCP_CSV, spec)
    assert [p.name for p in pts] == ["GCP_1", "2", "3"]
    assert pts[2].x == 500070.125 and pts[0].z == 12.345
    assert any("renamed" in w for w in warns)


def test_gcp_parse_whitespace_with_sigma():
    text = "# comment\nA 1 2 3 0.01\nB 4 5 6 0.02\n"
    spec = gcpmod.guess_spec(text)
    assert spec.delimiter is None and spec.skip_rows == 0
    spec.columns["sigma"] = 4
    pts, _ = gcpmod.parse_points(text, spec)
    assert pts[1].sigma == 0.02


def test_gcp_missing_role():
    with pytest.raises(ValueError):
        gcpmod.parse_points("a 1 2 3\n", gcpmod.TableSpec(columns={"name": 0, "x": 1}))


@pytest.mark.skipif(mmvii.find_mmvii() is None, reason="MMVII binary not available")
def test_import_with_mmvii_roundtrip(tmp_path):
    pts, _ = gcpmod.parse_points(synth.GCP_CSV, gcpmod.guess_spec(synth.GCP_CSV))
    gcpmod.import_with_mmvii(tmp_path, "Ctrl", pts, "EPSG:6653", sigma=0.02)
    sets = gcpmod.list_mmvii_sets(tmp_path)
    assert [s.name for s in sets] == ["Ctrl"]
    s = sets[0]
    assert s.crs == "EPSG:6653"
    assert [(g.name, g.x, g.z) for g in s.points] == [(p.name, p.x, p.z) for p in pts]
    assert s.points[0].sigma == pytest.approx(0.02)


@pytest.mark.skipif(mmvii.find_mmvii() is None, reason="MMVII binary not available")
def test_import_with_mmvii_converts_crs(tmp_path):
    pts = [gcpmod.GCP("P", -123.0, 49.25, 10.0)]
    gcpmod.import_with_mmvii(tmp_path, "Ctrl", pts, "EPSG:4326", "EPSG:32610")
    s = gcpmod.list_mmvii_sets(tmp_path)[0]
    assert s.crs == "EPSG:32610"
    assert s.points[0].x == pytest.approx(500000, abs=1)


def test_north_first():
    assert gcpmod.north_first("EPSG:4326")
    assert not gcpmod.north_first("EPSG:6653")
    assert not gcpmod.north_first("EPSG:32610")


def test_missing_grids_reports_geoid_or_nothing():
    grids = __import__("mmvii_studio.project", fromlist=["missing_grids"]).missing_grids(
        CRS.from_epsg(4979), CRS.from_epsg(6653))
    assert grids == [] or any("CGG2013" in g for g in grids)


def test_micasense_without_xmp_uses_file_name(tmp_path):
    """Early RedEdge firmware (e.g. sample.micasense.com example 1) writes no XMP, Make or Model."""
    from PIL import Image
    for idx in (1, 4):
        info = synth._tiff_info("", "", 36.67, -117.56, 1205.657, 4.2)
        del info[0x010F], info[0x0110]
        Image.new("I;16", (16, 12)).save(tmp_path / f"IMG_0455_{idx}.tif", tiffinfo=info)
    p = Project.open(tmp_path)
    p.scan_images()
    assert [(Path(im.path).name, im.band, im.band_index) for im in p.images] == [
        ("IMG_0455_1.tif", "Blue", 1), ("IMG_0455_4.tif", "NIR", 4)]
    assert p.images[0].make == "" and p.images[0].alt == pytest.approx(1205.657)


# --- MMVII camera database and metadata rules -------------------------------------

def test_fmt_writes_plain_decimals():
    # MMVII parses its command line itself and does not accept exponent notation, so
    # small pixel pitches must come out as plain decimals.
    assert mmvii._fmt(5) == "5"
    assert mmvii._fmt(22.0) == "22"
    assert mmvii._fmt(8.7606) == "8.7606"
    assert mmvii._fmt(0.0000015) == "0.000002"  # rounded, never "1.5e-06"
    assert "e" not in mmvii._fmt(1e-7)


def test_set_metadata_rule_rejects_an_unknown_mode(tmp_path):
    with pytest.raises(ValueError):
        mmvii.set_metadata_rule(tmp_path, "FocalMM", "22", ".*")


@pytest.mark.skipif(mmvii.find_mmvii() is None, reason="MMVII binary not available")
def test_register_camera_writes_the_local_database(tmp_path):
    mmvii.register_camera(tmp_path, "Studio Test Cam", (8.7606, 8.7606),
                          (22.812, 15.224), (2604, 1738))
    db = tmp_path / "MMVII-PhgrProj" / "MetaData" / "Std" / "CameraDataBase.xml"
    assert db.is_file(), "EditCamDataBase should write the project-local database"
    text = db.read_text()
    assert "Studio Test Cam" in text
    assert "2604" in text and "1738" in text


@pytest.mark.skipif(mmvii.find_mmvii() is None, reason="MMVII binary not available")
def test_register_plane_camera_writes_camera_and_rules(tmp_path):
    from mmvii_studio.raw import RawInfo, bayer_camera

    info = RawInfo(path="IMG_1753.CR2", width=5208, height=3476,
                   planes=["R", "G1", "G2", "B"], cfa_pattern="RGGB")
    cam = bayer_camera(info, 25.4 / 5798.66, "Canon EOS M", focal_mm=22.0)
    mmvii.register_plane_camera(tmp_path, cam)

    std = tmp_path / "MMVII-PhgrProj" / "MetaData" / "Std"
    assert cam.name in (std / "CameraDataBase.xml").read_text()
    rules = (std / "CalcMTD.xml").read_text()
    # Both the camera name and the focal length are attached to the plane files, because
    # MMVII reads no EXIF from TIFF at all and the planes would otherwise have no camera.
    assert f'<Val>"{cam.name}"</Val>' in rules
    assert '<Mode>"Focalmm"</Mode>' in rules and '<Val>"22"</Val>' in rules
    assert cam.file_pattern in rules


@pytest.mark.skipif(mmvii.find_mmvii() is None, reason="MMVII binary not available")
def test_registered_camera_lets_oricreatecalib_run(tmp_path):
    """The point of the camera database: OriCreateCalib on files with no readable EXIF."""
    from mmvii_studio.raw import RawInfo, bayer_camera

    info = RawInfo(path="IMG_1753.CR2", width=64, height=48,
                   planes=["R", "G1", "G2", "B"], cfa_pattern="RGGB")
    cam = bayer_camera(info, 25.4 / 5798.66, "Canon EOS M", focal_mm=22.0)
    for plane in ("R", "G1", "G2", "B"):
        Image.new("I;16", cam.image_size_px, 2048).save(tmp_path / f"IMG_1753_{plane}.tif")
    mmvii.register_plane_camera(tmp_path, cam)

    mmvii.run(["OriCreateCalib", r"IMG_1753_G1\.tif", "CalibG1"], cwd=tmp_path)
    calib = list((tmp_path / "MMVII-PhgrProj" / "Ori" / "CalibG1").glob("*.xml"))
    assert calib, "OriCreateCalib should write a calibration for the plane camera"
    # MMVII names the calibration after the camera and the focal length in microns.
    assert calib[0].name == "CalibIntr_CamCanon_EOS_M_Bayer_Foc22000.xml"
    text = calib[0].read_text()
    # The focal in pixels is the whole point of registering the plane geometry: MMVII
    # has to reach the same number bayer_camera does from the doubled pixel pitch.
    assert f"<NbPix>{cam.image_size_px[0]} {cam.image_size_px[1]}</NbPix>" in text
    f = float(re.search(r"<F>([\d.]+)</F>", text).group(1))
    # Not exact: _fmt gives MMVII the pixel pitch to six decimals, so MMVII divides by a
    # slightly rounded 8.760645 micron. That costs about 3e-8 relative, far inside any
    # real calibration, but it is why this is not an equality.
    assert f == pytest.approx(cam.focal_px, rel=1e-6)
    assert f == pytest.approx(2511.23, abs=0.01)


@pytest.mark.skipif(mmvii.find_mmvii() is None, reason="MMVII binary not available")
def test_every_metadata_mode_is_one_mmvii_accepts(tmp_path):
    """MTD_MODES is a whitelist, so a wrong spelling there turns a clear ValueError into
    an opaque MMVII failure. These are MMVII's serialised names, not its enum
    identifiers: eNbPixel is spelled "NbPix"."""
    for mode in mmvii.MTD_MODES:
        mmvii.set_metadata_rule(tmp_path, mode, "1", "XXXXXXXXXX.*XXXXXXX")
    modes = (tmp_path / "MMVII-PhgrProj" / "MetaData" / "Std" / "CalcMTD.xml").read_text()
    for mode in mmvii.MTD_MODES:
        assert f'<Mode>"{mode}"</Mode>' in modes
