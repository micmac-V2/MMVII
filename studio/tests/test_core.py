import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
import synth  # noqa: E402

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
