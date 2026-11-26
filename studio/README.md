# MMVII Studio

A desktop front end for MMVII, aimed at working through drone-mapping tutorials such as
E. Saczuk, *Processing Multi-spectral Imagery with Agisoft MetaShape Pro* (BCIT), with MMVII
doing the photogrammetry. Studio runs the `MMVII` binary and displays what it writes in
`MMVII-PhgrProj/`.

## Phase 1 (this version)

| Book step | Studio |
|---|---|
| Import photos | *Project ▸ Open folder* and *Import ▸ Scan images* read EXIF and XMP: GPS, camera, focal length, pixel size, band |
| Show only one band | The **Show** box filters by band. MicaSense `IMG_xxxx_N.tif` and DJI Mavic 3M `_MS_*.TIF` / `_D.JPG` files are grouped into captures |
| Camera centres (Fig. 2) | **Map** tab: camera centres and GCPs in the project CRS. Click a point to select its image |
| Reference pane | **Images** and **Ground control points** docks |
| Convert reference / set CRS (EPSG:6653) | *Project ▸ Coordinate system…* (empty = UTM zone of the GPS positions) |
| Import Reference (CSV, start at row 2) | *Import ▸ Ground control points…*: column roles, header rows, GCP CRS, accuracy. Runs `MMVII ImportGCP` and converts to the project CRS when needed |
| Geoids (ch. 1.1) | The **3D** tab warns when PROJ is missing the geoid grid needed for heights |
| Camera calibration info | **Cameras** tab: one row per camera model, band and image size |
| (not in the book) Camera raw | *Import ▸ Extract raw colour bands…* splits a CR2/NEF/ARW/DNG Bayer mosaic into one 16-bit image per colour channel, and registers a camera for them with MMVII |
| (not in the book) Camera motion | *Import ▸ BRAW camera motion…* reads the IMU out of Blackmagic RAW clips and writes an attitude per frame, plus the rolling-shutter readout time |

Not yet available: marker placement in images, alignment and bundle adjustment, camera
calibration estimation, dense cloud, DEM, orthomosaic.

## Run

Needs a built MMVII (https://github.com/micmac-V2/MMVII). From this folder (the repo root
of MMVII-studio, or `studio/` inside an MMVII checkout):

```bash
uv venv .venv && uv pip install --python .venv/bin/python -e '.[test,raw]'
.venv/bin/python -m mmvii_studio [project-folder]
```

The `raw` extra pulls in LibRaw (through rawpy) and is needed only to split raw files
into colour bands; everything else works without it.

Studio looks for MMVII in `$MMVII_BIN`, then on `PATH`, then in `../bin/MMVII` (when
`studio/` sits inside an MMVII checkout).

The project folder is the folder holding the images (subfolders such as MicaSense
`0000SET/000/` are scanned too). Studio writes `studio.json` and, through MMVII,
`MMVII-PhgrProj/` there.

## Camera raw and colour depth

MMVII computes on single-channel images: `cIm2D` with `eTN_U_INT2` (16-bit) or
`eTN_REAL4`. Its `cRGBImage` is three 8-bit planes and is marked in
`include/MMVII_Image2D.h` as "made for visu, not for intensive computation", so colour
is not something MMVII can match on. A raw file handed to MMVII directly is read as the
**8-bit JPEG preview** inside it, which throws away the sensor's full depth: a Canon EOS
M CR2 holds about 13 000 distinct levels per channel against 256 in that preview.

So Studio does not demosaic. *Import ▸ Extract raw colour bands…* splits the colour
filter array into one single-channel plane per filter — `IMG_1753_R.tif`, `_G1`, `_G2`,
`_B` — each half width and half height, 16-bit, which MMVII reads as `U_INT2`. The four
planes of one shot become the bands of one capture, exactly as a MicaSense rig's five
files do, and the raw file stays in the capture as the band `Raw`. Nothing is
interpolated: demosaicing invents detail that is correlated between channels and can
bias tie points, while a split mosaic keeps every value a real measurement. The two
greens are kept apart rather than averaged.

Each plane is written with the camera, the focal length and **half** the sensor's
focal-plane resolution, since its pixels cover two sensor pixels. Studio and exiftool
read that back; MMVII does not, because it reads no EXIF from TIFF files at all (a
synthetic TIFF whose EXIF exiftool reads fine still shows "No Exif metadata" in
`MMVII ImageMetadata`). MMVII therefore gets the geometry two other ways, both written
by *Extract raw colour bands…* into the project and not into the shared MMVII
installation:

- a camera in `MMVII-PhgrProj/MetaData/Std/CameraDataBase.xml` (the `LocalFolder`
  database, which MMVII loads last and so overrides `MMVII-RessourceDir/`)
- `ModelCam` and `Focalmm` rules in `MMVII-PhgrProj/MetaData/Std/CalcMTD.xml`, matched
  against the plane filenames

After that, `MMVII OriCreateCalib "IMG_...._G1.tif" CalibG1` succeeds, with the focal
length in pixels halved against the full frame (22 mm / 8.7606 µm = 2511 px for the
EOS M and its EF-M 22 mm).

## Blackmagic RAW camera motion

A Blackmagic Pocket Cinema Camera (firmware 7.9 and later) records a 6-DOF IMU into the
.braw file: a gyroscope in rad/s and an accelerometer in m/s². Neither needs Blackmagic's
SDK to read — a .braw is an ISO-BMFF container, the clip metadata is a QuickTime
`keys`/`ilst` table in `moov/meta`, and the IMU lives in two timed-metadata tracks whose
samples are 20-byte boxes (`mogy` for gyroscope, `moac` for accelerometer, three
little-endian float32 each).

*Import ▸ BRAW camera motion…* reads every .braw in the project and writes, next to each
clip, a `<clip>_attitude.csv` (per frame: roll, pitch, yaw and a quaternion) and a
`<clip>.ori` for `MMVII ImportOri`:

```bash
MMVII ImportOri A004_01010039_C042.ori NXYZWPK Calib InitIMU AngU=degree
```

Two things it gives you that the pictures alone do not:

- **Attitude per frame.** The accelerometer sees gravity, so roll and pitch do not drift.
  Yaw is *not* observable from a gyroscope and accelerometer alone: it starts at zero and
  drifts, and fixing heading needs a magnetometer, GNSS motion, or the photogrammetry.
- **Rolling-shutter timing.** `sensor_line_time` x the captured rows gives the readout
  time — 7.83 ms for the Pocket 4K's 1920x1080 crop, 16.22 ms at 4096x2160. On a moving
  CMOS camera that is a systematic geometric error in every frame.

Per-frame `bmdf` boxes also carry focal length (`"14mm"`), aperture, ISO, shutter angle,
focus distance and white balance. Focal length matters most: on a zoom it can change
during a clip. The sensor *size* is not in the file, so `KNOWN_SENSORS` in
[braw.py](mmvii_studio/braw.py) holds published figures (Pocket 4K: 18.96 mm across 4096
native pixels) to turn mm into pixels — 14 mm / 4.6289 µm = 3024 px. Verify those against
your own body before trusting a calibration built on them.

### Things found the hard way

- **Firmware 7.9.1 lies about the IMU sample rate.** The metadata tracks always declare a
  timescale of 999 with a sample delta of 1, i.e. 999 Hz, whatever the real rate is.
  Measured against the video duration on one Pocket 4K: 999 Hz at 24.000 fps, but 395.6 Hz
  at 23.976 fps, 664.1 Hz at 29.97 fps and 1030.9 Hz at 60 fps. Firmware 8.1 clips were
  correct. Believing the declared rate on a 23.976 fps clip stretches IMU time by 2.5x and
  ruins the integration, so `read_imu` derives the rate from the video duration whenever
  the two disagree by more than 2%, and says so in the log.
- **OIS does not have to be off.** Clips recorded with `ois_enable = 1` still carried a
  full set of IMU samples, contrary to a widely repeated forum claim.
- **The accelerometer reads about 2.4% high.** Gravity came out at 10.01–10.07 m/s² across
  clips instead of 9.807, consistently, which is a scale error rather than noise. It does
  not affect the gravity *direction*, so attitude is unharmed, but do not use the
  magnitudes without calibrating.
- **The positions in the .ori file are zeros.** An IMU gives attitude, not position;
  integrating its accelerometer twice drifts hopelessly. MMVII's format demands X, Y and Z,
  so they are placeholders to be replaced from GNSS or resection, or left for the bundle
  adjustment. The angles likewise need a boresight rotation between IMU and lens before
  they are camera angles — on the clips tested, gravity lies along the IMU's Y axis, which
  shows up as a constant ~90° roll if you ask for `reference="gravity"`.

Decoding pictures out of a .braw still needs Blackmagic's SDK, so the Bayer-plane
splitting described above does not apply to BRAW.

## Notes

- **Axis order.** MMVII reads coordinates in the CRS's official axis order, which is
  latitude first for EPSG:4326. Studio keeps x = easting/longitude and y = northing/latitude,
  and swaps them when writing files for MMVII.
- **GCP files.** MMVII splits fields on spaces and commas, so empty CSV fields would shift
  columns. Studio writes a clean `GCP-<set>.txt` file and imports that instead.

## Tests

```bash
.venv/bin/python -m pytest tests
.venv/bin/python tests/make_demo.py /tmp/demo   # synthetic RedEdge-M flight + Boundary.csv
```

## Sample data

- MicaSense RedEdge raw captures (about 8 MB each): https://sample.micasense.com
- ICAERUS DJI Mavic 3M and Mavic 3 datasets: https://github.com/ICAERUS-EU/Zenodo_Datasets
  (for example https://zenodo.org/records/15114347, 342 RGB images, 6.1 GB)

## Licence

CeCILL-B, the same licence as MMVII. See [LICENSE.md](LICENSE.md).
