# Handoff: Metashape reference runs on Windows

Written 2026-10-14 at the end of a Claude Code session on the Mac (`~/code/MMVII`). Give
this file to the person or Claude session on the Windows machine. It assumes no other context.

## The project in one paragraph

**MMVII Studio** is a PySide6 + pyvista desktop front end for the MMVII photogrammetry
tools (https://github.com/micmac-V2/MMVII). The aim is to work through E. Saczuk,
*Processing Multi-spectral Imagery with Agisoft MetaShape Pro* (BCIT; MicaSense RedEdge-M,
5 bands, GCPs in EPSG:6653) with MMVII instead of Metashape. Agisoft Metashape Pro is the
**reference**: we run the same data through it and compare its calibration, poses and GCP
errors with MMVII's. Studio runs the `MMVII` binary and displays what it writes in
`MMVII-PhgrProj/`.

## Where the code is

| Place | What |
|---|---|
| https://github.com/KevinGriffin-new/MMVII-studio (private), `main` | Studio only (the `studio/` folder of the MMVII branch, history kept) |
| `kevin@fedora.local:git/MMVII.git`, branch `feat/mmvii-studio` | Full MMVII fork with `studio/` inside |
| Mac: `~/code/MMVII`, branch `feat/mmvii-studio` | Working copy. Studio commits go to GitHub with `git subtree push --prefix=studio studio main` |

## Status

**Phase 1 is done** (commit `afeb0baf` in MMVII-studio):
- scans image EXIF/XMP and groups MicaSense `IMG_xxxx_N.tif` and DJI Mavic 3M `_MS_*` files into captures and bands
- imports GCPs through `MMVII ImportGCP`
- map, 3D and camera-summary views
- 14 tests

It was tested on synthetic flights and on MicaSense sample #1 (https://sample.micasense.com,
which is a single capture).

**Phase 2, camera calibration, is next.** It needs:
- a real flight with overlap
- tie points
- the Metashape reference, which is this handoff's job

**Since then**, Blackmagic RAW camera-motion ingest also landed: `mmvii_studio/braw.py`
reads a .braw container directly (no Blackmagic SDK) for clip metadata, the gyroscope and
accelerometer tracks, and per-frame focal length, and writes a per-frame attitude file that
`MMVII ImportOri` accepts. Verified end to end on a Pocket Cinema Camera 4K: MMVII's own
reconstructed quaternion matched ours to ~1e-6. Watch out for firmware 7.9.1, which
misdeclares the IMU sample rate; the README's "Things found the hard way" has the numbers.

**Also since then**, camera raw ingest landed on the Mac (see the "Camera raw and colour depth"
section of `studio/README.md`): CR2 and friends are scanned through a TIFF-IFD reader, a Bayer
mosaic can be split into four 16-bit single-channel band images, and the camera-name problem
below is solved for them. The test camera is a Canon EOS M with an EF-M 22 mm, which has **no
GPS**, so it gives calibration and tie-point material but no camera positions.

## Facts found so far (don't re-derive)

- **Demo mode can't produce the reference.** On the Mac, Metashape Pro 2.3.2 runs in demo
  mode. A probe script found that scripting, adding photos and processing work, but
  `chunk.exportCameras()` and `doc.save()` fail with "No license found". Use the licensed
  Windows copy, registered with the BCIT licence server `agisoftgis.csg.bcit.ca:5842` over
  the school VPN.
- **MMVII uses the CRS's official axis order.** For EPSG:4326 a GCP file must be lat, lon.
  Studio keeps x = easting/longitude and swaps when needed. Projected systems such as
  EPSG:6653 and EPSG:326xx are easting first.
- **Geoid grid.** PROJ on the Mac lacks `ca_nrc_CGG2013i08.tif` (CGVD2013). Without it,
  EPSG:6653 heights fall back to ellipsoidal heights with no error. This matters for comparing
  heights with Metashape, which ships its own geoids folder. Install it with
  `projsync --file ca_nrc_CGG2013i08.tif`.
- **MMVII needs a camera name, and cannot read EXIF from TIFF.** `MMVII OriCreateCalib` fails
  with "Camera Name is not init" on TIFF input. The reason recorded earlier — that the RedEdge
  files have no Make or Model — is not the whole story: **MMVII reads no EXIF from TIFF files
  at all.** A synthetic TIFF carrying Make, Model and FocalLength that exiftool reads correctly
  still gives "No Exif metadata" from `MMVII ImageMetadata`, while the same tags in a JPEG are
  read. So any TIFF workflow (MicaSense, DJI MS, split Bayer planes) must supply the camera
  explicitly, whatever the files contain. MMVII's `CameraDataBase.xml` also has no MicaSense or
  DJI entries.

  This is now done for camera raw, and the same two steps work for MicaSense and DJI
  (`mmvii_studio/mmvii.py`):
  - `EditCamDataBase <name> LocalFolder 1 Param=[SzPixX,SzPixY,SzSensX,SzSensY,SzImX,SzImY]`
    writes the geometry to `MMVII-PhgrProj/MetaData/Std/CameraDataBase.xml`. The enum value is
    `LocalFolder`, not `eLocalFolder`. MMVII loads that database **last**, so it overrides the
    shared `MMVII-RessourceDir/CameraDataBase.xml` without touching that tracked file.
  - `EditCalcMTDI ModelCam Value=<name> PatternApply=<regex> Save=true`, and the same for
    `Focalmm`, write per-project rules to `MMVII-PhgrProj/MetaData/Std/CalcMTD.xml`.
    The mode names are the strings in the `eMTDIm` table of `src/Serial/ElemStrToVal.cpp`,
    not the enum identifiers in `include/MMVII_enums.h` with the `e` dropped: all match
    except `eNbPixel`, which MMVII spells **`NbPix`** on the command line.

  Verified end to end on Canon EOS M CR2 files: `OriCreateCalib` then builds
  `CalibIntr_CamCanon_EOS_M_Bayer_Foc22000` with F = 2511.23 px, PP at the image centre.
- **MMVII has no tie-point generator of its own in normal use.** Its tutorial gets tie points
  from MicMac v1 (`mm3d Tapioca`), which isn't installed. Undecided: build MicMac v1, use
  MMVII's experimental `TieP-Aime*`, or import Metashape tie points (`MMVII ImportTiePMul`).

## Task on Windows

### 0. Prerequisites
- School VPN up. Metashape Pro licensed: Help ▸ Activate Product shows the floating licence,
  or check that `Metashape.app.activated` is `True` in the Python console.
- A dataset with overlapping images. Candidates, not chosen yet; confirm the size before
  downloading:
  - a RedEdge flight of your own with surveyed GCPs (best match for the book)
  - an ICAERUS DJI Mavic 3M flight: https://github.com/ICAERUS-EU/Zenodo_Datasets
    (Canyelles vineyard records, e.g. https://zenodo.org/records/15064455). These have no GCPs.

### 1. Write `reference/metashape_reference.py` (in MMVII-studio)
A headless script, run as:
```
"C:\Program Files\Agisoft\Metashape Pro\metashape.exe" -r reference\metashape_reference.py <images_dir> <out_dir> [gcp.csv] [EPSG]
```
It should:
1. Stop with a clear message if `Metashape.app.activated` is false.
2. For each band separately, as the book does (filter by `IMG_xxxx_<band>.tif` or
   `_MS_<band>`), make one chunk, add the photos and record the Metashape version and
   parameters.
3. Follow the book's chapter 1.2 settings, which were checked against Fig. 6:
   - Estimate Image Quality and leave out photos scoring below 0.6. (In 2.x this is probably
     `chunk.analyzeImages()` and the `Image/Quality` metadata key; check.)
   - Convert the camera reference to the project CRS (the book uses EPSG:6653).
   - Align Photos with:
     - Accuracy **Highest** (`matchPhotos(downscale=0)`)
     - Generic preselection **off**
     - Reference preselection **on**, mode **Source**
     - Key point limit **40 000**, tie point limit **4 000**
     - Guided image matching **off**, masks none
     - Adaptive camera model fitting **off** (`alignCameras(adaptive_fitting=False)`)
   - Then `optimizeCameras` with the default parameters (the book's "Optimize Cameras, accept
     the defaults").
4. If a GCP CSV is given: set `chunk.crs` (e.g. EPSG:6653) and import the markers.
   Projections can't be placed headless, so only export reference positions and, if they
   exist, Metashape's marker errors.
5. Export into `<out_dir>/<band>/`:
   - `cameras.xml`: `chunk.exportCameras(..., format=Metashape.CamerasFormatXML)` (poses and calibration)
   - `calibration.xml`: `sensor.calibration.save(..., format=Metashape.CalibrationFormatXML)` for each sensor
   - `reference.csv`: camera and marker estimated positions and errors (`chunk.exportReference`)
   - `report.pdf`: `chunk.exportReport`
   - tie points in a format MMVII can import, if one exists (look at `exportCameras` formats
     such as Bundler, or `exportPointCloud` with tie-point source). Write down what was possible.
   - `run.json`: Metashape version, all parameters, image count, aligned count, timings
6. Optionally, a second run with `addPhotos(..., layout=Metashape.MultiplaneLayout)` (all
   bands as one rig) to get each band's calibration plus the band-to-band offsets.

Check every API name against the Metashape 2.x Python reference for the installed version.
The names above are from memory.

### 2. Bring the results back
Commit the output folder to MMVII-studio on a branch `reference/<dataset>`. The XML, CSV and
JSON files are small, but leave out images and `.psx` projects. Note in `run.json` which
dataset and version were used.

### 3. Back on the Mac (not part of this handoff)
- a Metashape-calibration → MMVII converter (f, cx, cy, k1–k4, p1, p2, b1, b2 in pixels, with
  Metashape's principal point measured from the image centre)
- Studio's calibration view comparing Metashape and MMVII
- MMVII `OriCreateCalib` and `OriBundleAdj` on the same images

## Conventions
- Keep commits small. Commit messages end with a `Co-Authored-By` line when Claude wrote the commit.
- Reference outputs never go into the MMVII fork, only into MMVII-studio.
- Studio code is CeCILL-B (same as MMVII).
