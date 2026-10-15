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

Not yet available: marker placement in images, alignment and bundle adjustment, camera
calibration estimation, dense cloud, DEM, orthomosaic.

## Run

```bash
cd studio
uv venv .venv && uv pip install --python .venv/bin/python -e '.[test]'
.venv/bin/python -m mmvii_studio [project-folder]
```

Studio looks for MMVII in `$MMVII_BIN`, then on `PATH`, then in `../bin/MMVII`.

The project folder is the folder holding the images (subfolders such as MicaSense
`0000SET/000/` are scanned too). Studio writes `studio.json` and, through MMVII,
`MMVII-PhgrProj/` there.

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
