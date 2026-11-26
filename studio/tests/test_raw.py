"""Tests for camera-raw ingest: CFA splitting, raw EXIF reading, plane camera geometry."""

import sys
from pathlib import Path

import numpy as np
import pytest
from PIL import Image

sys.path.insert(0, str(Path(__file__).parent))
import synth  # noqa: E402

from mmvii_studio import raw as rawmod
from mmvii_studio.bands import classify, group_captures, bands_of
from mmvii_studio.exif import IMAGE_SUFFIXES, ImageInfo, read_image_info, read_tiff_tags

RGGB = np.array([[0, 1], [3, 2]])     # as LibRaw reports a Canon sensor
RGGB_DESC = b"RGBG"                   # index 0=R, 1=G, 2=B, 3=the second G


def test_is_raw_and_suffixes():
    assert rawmod.is_raw("IMG_1753.CR2") and rawmod.is_raw(Path("a/b/x.nef"))
    assert not rawmod.is_raw("IMG_1753.JPG") and not rawmod.is_raw("plane.tif")
    # raw files must be picked up by the scanner alongside rendered ones
    assert ".cr2" in IMAGE_SUFFIXES and ".jpg" in IMAGE_SUFFIXES


def test_plane_names_disambiguates_the_two_greens():
    pattern, names = rawmod._plane_names(RGGB, RGGB_DESC)
    assert pattern == "RGGB"
    assert names == ["R", "G1", "G2", "B"]


def test_split_mosaic_strides_each_plane_correctly():
    m = synth.bayer_mosaic(width=8, height=6, base=1000)
    planes = rawmod.split_mosaic(m, RGGB, RGGB_DESC, black_levels=[0, 0, 0, 0],
                                white_level=65535, normalize=False)
    assert sorted(planes) == ["B", "G1", "G2", "R"]
    for p in planes.values():
        assert p.shape == (3, 4) and p.dtype == np.uint16
    # synth gives plane (dy,dx) a base of 1000*(1+dy*2+dx); R is (0,0), G1 (0,1),
    # G2 (1,0), B (1,1). A wrong stride would pick up a different base.
    assert planes["R"][0, 0] == 1000
    assert planes["G1"][0, 0] == 2000
    assert planes["G2"][0, 0] == 3000
    assert planes["B"][0, 0] == 4000


def test_split_mosaic_normalize_maps_black_to_zero_and_white_to_full():
    m = np.array([[2047, 16383], [2047, 16383]], dtype=np.uint16)
    planes = rawmod.split_mosaic(m, RGGB, RGGB_DESC, black_levels=[2047] * 4,
                                white_level=16383, normalize=True)
    assert planes["R"][0, 0] == 0          # at black level
    assert planes["G1"][0, 0] == 65535     # at white level
    # unnormalized keeps the sensor ADU untouched
    adu = rawmod.split_mosaic(m, RGGB, RGGB_DESC, black_levels=[2047] * 4,
                              white_level=16383, normalize=False)
    assert adu["R"][0, 0] == 2047 and adu["G1"][0, 0] == 16383


def test_split_mosaic_clips_below_black():
    m = np.full((2, 2), 1000, dtype=np.uint16)  # below the 2047 black level
    planes = rawmod.split_mosaic(m, RGGB, RGGB_DESC, black_levels=[2047] * 4,
                                white_level=16383, normalize=True)
    assert planes["R"][0, 0] == 0  # clipped, not wrapped around


def test_split_mosaic_trims_odd_dimensions():
    m = synth.bayer_mosaic(width=9, height=7)
    planes = rawmod.split_mosaic(m, RGGB, RGGB_DESC, black_levels=[0] * 4,
                                white_level=65535, normalize=False)
    # 9x7 trims to 8x6, so every plane is exactly half of that
    for p in planes.values():
        assert p.shape == (3, 4)


def test_split_mosaic_subset_of_planes():
    m = synth.bayer_mosaic()
    planes = rawmod.split_mosaic(m, RGGB, RGGB_DESC, black_levels=[0] * 4,
                                white_level=65535, normalize=False, planes=["G1", "G2"])
    assert sorted(planes) == ["G1", "G2"]


def test_split_mosaic_rejects_bad_input():
    m = synth.bayer_mosaic()
    with pytest.raises(ValueError):
        rawmod.split_mosaic(m[np.newaxis], RGGB, RGGB_DESC, black_levels=[0] * 4,
                            white_level=65535)
    with pytest.raises(ValueError):  # not a 2x2 CFA
        rawmod.split_mosaic(m, np.zeros((3, 3), int), RGGB_DESC, black_levels=[0] * 4,
                            white_level=65535)
    with pytest.raises(ValueError):  # white below black
        rawmod.split_mosaic(m, RGGB, RGGB_DESC, black_levels=[5000] * 4, white_level=100)


def test_plane_file_pattern_matches_only_plane_files():
    import re
    pat = rawmod.plane_file_pattern(["R", "G1", "G2", "B"])
    assert re.fullmatch(pat, "IMG_1753_G1.tif")
    assert re.fullmatch(pat, "IMG_1753_R.tif")
    assert not re.fullmatch(pat, "IMG_1753.CR2")
    assert not re.fullmatch(pat, "IMG_1753_G1.jpg")


def test_bayer_camera_halves_resolution_and_focal_in_pixels():
    info = rawmod.RawInfo(path="x.CR2", width=5208, height=3476,
                          planes=["R", "G1", "G2", "B"], cfa_pattern="RGGB")
    pitch_mm = 25.4 / 5798.66          # Canon EOS M focal-plane resolution
    cam = rawmod.bayer_camera(info, pitch_mm, "Canon EOS M", focal_mm=22.0)

    assert cam.name == "Canon EOS M Bayer"
    assert cam.image_size_px == (2604, 1738)
    # one photosite per 2x2 block, so the pitch doubles
    assert cam.pixel_size_micron[0] == pytest.approx(pitch_mm * 2000.0)
    # ... over the same physical sensor area
    assert cam.sensor_size_mm[0] == pytest.approx(5208 * pitch_mm, rel=1e-9)
    assert cam.sensor_size_mm[1] == pytest.approx(3476 * pitch_mm, rel=1e-9)
    # focal in mm is unchanged, so focal in pixels halves
    assert cam.focal_px == pytest.approx(22.0 / (pitch_mm * 2), rel=1e-9)
    assert cam.focal_px == pytest.approx(2511.23, abs=0.01)


def test_bayer_camera_rejects_bad_pitch():
    info = rawmod.RawInfo(path="x.CR2", width=100, height=100, planes=["R"])
    with pytest.raises(ValueError):
        rawmod.bayer_camera(info, 0.0, "Cam")


def test_plane_pixel_size_mm():
    assert rawmod.plane_pixel_size_mm(0.0043803) == pytest.approx(0.0087606)
    assert rawmod.plane_pixel_size_mm(None) is None


def test_read_tiff_tags_reads_a_canon_like_container(tmp_path):
    path = synth.make_canon_like_tiff(tmp_path)
    ifd0, exif, gps = read_tiff_tags(path)
    assert ifd0[0x010F] == "Canon"
    assert ifd0[0x0110] == "Canon EOS M"
    assert exif[0x9003] == "2026:11:10 15:01:16"
    num, den = exif[0x920A]
    assert num / den == pytest.approx(22.0)


def test_read_tiff_tags_ignores_a_non_tiff(tmp_path):
    junk = tmp_path / "not.cr2"
    junk.write_bytes(b"this is not a TIFF container at all")
    assert read_tiff_tags(junk) == ({}, {}, {})


def test_raw_image_info_without_libraw(tmp_path, monkeypatch):
    """EXIF must still be reported when LibRaw cannot open the file."""
    path = synth.make_canon_like_tiff(tmp_path, name="IMG_9002.CR2")

    def boom(_):
        raise RuntimeError("no LibRaw here")

    monkeypatch.setattr(rawmod, "read_raw_info", boom)
    info = read_image_info(path, tmp_path)
    assert info.is_raw is True
    assert info.make == "Canon" and info.model == "Canon EOS M"
    assert info.focal_mm == pytest.approx(22.0)
    assert info.pixel_size_mm == pytest.approx(25.4 / 5798.66, rel=1e-4)
    assert info.bit_depth is None          # unknown without LibRaw
    assert info.width == 0                 # geometry unknown too


def test_cfa_planes_group_into_one_capture():
    ims = [ImageInfo(path=f"IMG_1753_{p}.tif") for p in ("R", "G1", "G2", "B")]
    caps = group_captures(ims)
    assert list(caps) == ["IMG_1753"]
    assert sorted(caps["IMG_1753"]) == ["Blue", "Green1", "Green2", "Red"]
    assert bands_of(ims) == ["Blue", "Green1", "Green2", "Red"]


def test_cfa_pattern_does_not_capture_other_names():
    # lower case, and names that merely end in a letter, stay single-band RGB
    assert classify(ImageInfo(path="some_b.tif")) == ("some_b", "RGB")
    assert classify(ImageInfo(path="holiday.jpg")) == ("holiday", "RGB")
    # MicaSense and DJI keep their own grouping
    assert classify(ImageInfo(path="IMG_0012_1.tif")) == ("IMG_0012", "Blue")
    assert classify(ImageInfo(path="DJI_20230609120000_0012_MS_NIR.TIF"))[1] == "NIR"


def test_plane_exif_halves_resolution_and_suffixes_the_model(tmp_path):
    src = synth.make_canon_like_tiff(tmp_path, xres_per_inch=5798.66, focal=22.0)
    info = rawmod.plane_exif(src)
    assert info[0x010F] == "Canon"
    assert info[0x0110] == "Canon EOS M Bayer"
    sub = info[0x8769]
    assert float(sub[0x920A]) == pytest.approx(22.0)       # focal in mm unchanged
    assert float(sub[0xA20E]) == pytest.approx(5798.66 / 2, rel=1e-4)  # pixels twice as wide
    assert sub[0xA210] == 2                                # still per inch


def test_plane_exif_custom_suffix(tmp_path):
    src = synth.make_canon_like_tiff(tmp_path)
    assert rawmod.plane_exif(src, model_suffix=" CFA")[0x0110] == "Canon EOS M CFA"


def test_plane_exif_roundtrips_through_a_written_plane(tmp_path):
    """A written plane must report half the sensor's pitch when read back."""
    src = synth.make_canon_like_tiff(tmp_path)
    plane = np.full((4, 6), 1234, dtype=np.uint16)
    dest = tmp_path / "IMG_9001_G1.tif"
    rawmod._save_u16_tiff(plane, dest, rawmod.plane_exif(src))

    got = read_image_info(dest, tmp_path)
    assert got.model == "Canon EOS M Bayer"
    assert got.focal_mm == pytest.approx(22.0)
    # sensor pitch is 25.4/5798.66; the plane's is twice that
    assert got.pixel_size_mm == pytest.approx(2 * 25.4 / 5798.66, rel=1e-4)
    assert np.array_equal(np.array(Image.open(dest)), plane)  # pixels untouched by the tags


def test_plane_exif_survives_a_file_without_exif(tmp_path):
    bare = tmp_path / "bare.cr2"
    bare.write_bytes(b"not a tiff")
    info = rawmod.plane_exif(bare)
    assert 0x010F not in info and 0x0110 not in info  # nothing invented


def test_raw_file_is_its_own_band_not_rgb():
    """The raw source is the capture's source, not one of its colour bands."""
    raw = ImageInfo(path="IMG_1753.CR2", is_raw=True)
    assert classify(raw) == ("IMG_1753", "Raw")
    # and it groups with the planes extracted from it
    ims = [raw] + [ImageInfo(path=f"IMG_1753_{p}.tif") for p in ("R", "G1", "G2", "B")]
    caps = group_captures(ims)
    assert list(caps) == ["IMG_1753"]
    assert sorted(caps["IMG_1753"]) == ["Blue", "Green1", "Green2", "Raw", "Red"]
