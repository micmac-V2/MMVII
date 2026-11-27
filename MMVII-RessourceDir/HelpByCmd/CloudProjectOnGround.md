# `CloudProjectOnGround` — Programmer Documentation

**Command name:** `CloudProjectOnGround`
**Source file:** [`src/Mesh/cCloudProjectOnGround.cpp`](../../src/Mesh/cCloudProjectOnGround.cpp)
**Main class:** `cAppliCloudProjectOnGround`
**Spec object:** `TheSpecCloudProjectOnGround` — *"Project of Pax/Depth maps to ground"*, category `eApF::Cloud`, input/output type `eApDT::Image`.

---

## Purpose

`CloudProjectOnGround` converts a per-pixel 2D cloud — either a **disparity map** (epipolar
geometry) or a **depth map** (image + depth geometry) — into a **rasterized DEM**
("bascule" of the cloud onto the ground). For each tile of the input it:

1. projects every valid pixel to a 3D ground point,
2. builds an implicit triangulated mesh over the pixel grid,
3. runs a Z-buffer to classify each triangle (visible / hidden / distorted / …),
4. rasterizes the visible triangles at a given GSD into a georeferenced elevation tile,
5. finally merges all tiles into one global DEM + Mask (+ optional Correlation) raster.

It is the middle step of the dense-matching pipeline:

```
EpipGenDenseMatch                CloudProjectOnGround              CloudFuser
(src/DenseMatch/EpipGenDenseMatch.cpp)  (this command)             (src/Mesh/cCloudFuser.cpp)
   writes DMEpipPairResult.xml  -->  reads it, writes DEM tile  -->  reads cCloudRaster xml,
                                     rasters + cCloudRaster xml      fuses per-pair DEMs
```

---

## Class architecture

```cpp
class cAppliCloudProjectOnGround : public cMMVII_Appli,
                                   public cAppliParseBoxIm<tREAL4>
```

* **`cMMVII_Appli`** — standard MMVII application plumbing (argument parsing,
  `cPhotogrammetricProject`, timing segments, parallel sub-process recall).
* **`cAppliParseBoxIm<tREAL4>`** (declared in
  [`include/MMVII_Image2D.h`](../../include/MMVII_Image2D.h)) — tiled parsing of a large
  image. The constructor configures tiles of **2000×2000** pixels with a **50×50**
  overlap and `ParalTiles=true`, so each tile is processed by a *parallel recall* of the
  command itself (`APBI_ExecAll()` forks one sub-process per tile; each sub-process
  reenters `Exe()` with `InsideParalRecall() == true`).

Supporting types:

| Type | Header | Role |
|---|---|---|
| `cDMEpipPairResult` | [`include/MMVII_CloudRaster.h`](../../include/MMVII_CloudRaster.h) | serialized *input* descriptor (written by `EpipGenDenseMatch`) |
| `cCloudRaster` | [`include/MMVII_CloudRaster.h`](../../include/MMVII_CloudRaster.h) | serialized *output* descriptor (consumed by `CloudFuser`) |
| `cZBuffer`, `eZBufRes`, `cResModeSurfD` | [`include/MMVII_ZBuffer.h`](../../include/MMVII_ZBuffer.h), [`src/Geom3D/cZBuffer.cpp`](../../src/Geom3D/cZBuffer.cpp) | triangle visibility computation |
| `cTriangulation3D<tREAL8>`, `cMeshTri3DIterator` | `MMVII_Triangles.h` / `MMVII_ZBuffer.h` | explicit mesh + iterator adapter for the Z-buffer |
| `cSensorImage` / `cSensorCamPC` | `MMVII_PCSens.h` | camera models (generic sensor vs. central perspective) |

---

## Command-line interface

### Mandatory arguments (`ArgObl`)

| # | Member | Description |
|---|---|---|
| 1 | `mNameDMEpipPairResult` | XML file (tagged, `cDMEpipPairResult`) describing the dense-matching result of one epipolar pair |
| 2 | `mPhProj.DPOrient().ArgDirInMand()` | input orientation directory in the MMVII specif (`Ori`) |
| 3 | `mPhProj.DPMeshDev().ArgDirOutMand()` | output *mesh-dev* directory, used for intermediate per-tile rasters |

### Optional arguments (`ArgOpt`, plus the shared `APBI_ArgOpt` tiling options)

| Name | Member | Default | Meaning |
|---|---|---|---|
| `NameResult` | `mNameResult` | `Dem_<Im1>.tif` | prefix/name for output files |
| `GSD` | `mGSD` | auto (`ComputeGSD()`) | ground sampling distance of the output raster |
| `DirOut` | `mOutDir` | `VISU/CloudProjectOnGround/` | directory of the final merged outputs |
| `MII` | `mMII` | `0.0` | *Margin Inside Image* passed to `cSetVisibility` for triangle validation |
| `OriInFatherDir` | `mOriInFatherDir` | `false` | look up the orientation directory in the *father* project directory instead of `DirProject()` (see §5.1) |
| `DistortionThreshold` | `mStretschingThresh` | `4.0` | intended max triangle stretching before a triangle is discarded |

`APBI_ArgOpt` additionally exposes the standard tiling options of `cAppliParseBoxIm`
(test box, tile size, overlap…).

### Example (from `Samples()`)

```
MMVII CloudProjectOnGround /path/to/DMEpipPairResult.xml ori BascTmpFiles GSD=0.1
```

### Input XML structure (`cDMEpipPairResult`)

```xml
<Data>
    <ModeGeom>DEPTH</ModeGeom>              <!-- "EPIP" or "DEPTH" (eModeGeom) -->
    <NameEpipIm1>Epip1.png</NameEpipIm1>
    <NameEpipIm2>Epip2.png</NameEpipIm2>
    <NameDisparity>DenseMatchResult/disparity.tif</NameDisparity>
    <NameCorrel>DenseMatchResult/Correlation.tif</NameCorrel>   <!-- optional -->
    <NameMasq>DenseMatchResult/occlusion_mask.tif</NameMasq>
</Data>
```

All image paths are relative to `DirProject()`. If `NameCorrel` is non-empty,
`mBascCorrel` is set and a correlation raster is also produced.

---

## Geometry modes (`eModeGeom`)

Defined in [`include/MMVII_enums.h`](../../include/MMVII_enums.h):

* **`eGEOM_EPIP`** — the 2D cloud is a *disparity map* between two epipolar images.
  A pixel `(x, y)` of image 1 with disparity `d` matches `(x+d, y)` in image 2; the 3D
  point is recovered with `RobustBundleInters()` on the two `Image2Bundle()` rays. Both
  sensors (`mCamGen`, `mSecCamGen`) are required.
* **`eGEOM_DEPTH`** — the 2D cloud is a *depth map*; the 3D point is
  `mCamGen->ImageAndDepth2Ground()`. A depth value of **0 denotes no-data** and the pixel
  is skipped. The camera must be a `cSensorCamPC` (asserted).

---

## Execution flow

### `Exe()` — master and recall entry point

```
Exe()
 |- (OriInFatherDir) DPOrient().SetAllowDirInEmpty()
 |- mPhProj.FinishInit()
 |- (OriInFatherDir) redirect DPOrient full input dir to <father>/MMVII-PhgrProj/Ori/<OriIn>/
 |- ReadFromFile_Std(cDMEpipPairResult, mNameDMEpipPairResult)   -> mode, image/masq/correl names
 |- resolve mOutDir (default DirVisuAppli()) and mNameResult (default "Dem_<Im1>.tif")
 |- mCamGen = mPhProj.ReadSensor(Im1);  (EPIP) mSecCamGen = ReadSensor(Im2)
 |- if GSD not user-set:  mGSD = ComputeGSD(); push "GSD=<val>" into mArgsRecall
 |        (so every parallel tile recall reuses the master's estimate instead of recomputing)
 |- APBI_ExecAll()               -> for each tile: ExeOnParsedBox()  [in parallel recalls]
 `- if !InsideParalRecall():     -> master only
      |- MergeResults()
      |- RemovePatternFile("BLOC.*" / "MASQ.*" / "CORR.*")   (delete tile intermediates)
      `- SaveInFile(cCloudRaster, <Prof_...>.xml)            (descriptor for CloudFuser)
```

**`OriInFatherDir`:** the epipolar project usually lives in a sub-folder of the real
project; with this flag the `Ori` directory is resolved against the directory *above*
`MMVII-PhgrProj` instead of `DirProject()`.

### `ExeOnParsedBox()` — per tile

Called once per tile by `APBI_ExecAll()`:

1. `mImPx1 = APBI_ReadIm<tREAL4>(disparity/depth)` — tile crop of the 2D cloud,
2. `mImMasq1 = ReadMasqWithDef(CurBoxIn(), masq)` — validity mask (all-1 if the file is missing),
3. if `NameCorrel` is set: read `mImCorrel`, set `mBascCorrel = true`,
4. `MakeFastBasc()`.

### `BascOnePoint(A, offset, &valid)` — pixel -> 3D

Core geometric routine; returns the pair
`{aPWithDepth, aPZA}` = `{(x, y, depth-like) in image space, 3D ground point}`.

* **EPIP:** pseudo-intersection of the two epipolar rays. The "depth-like" value stored
  with the 2D point is `-Z` for a generic sensor, or the true
  `Ground2ImageAndDepth(...).z()` for a `cSensorCamPC`.
* **DEPTH:** direct `ImageAndDepth2Ground`; `depth == 0` sets `valid=false` (no-data).

Side effect: adds `Proj(aPZA)` to `mBoxLocTarget`, the ground footprint accumulator for
the current tile (used later to size and georeference the output tile raster).

### `ComputeGSD()` — automatic GSD estimation

Only runs when the user didn't pass `GSD=`. It reads a small **100×100** window at the
center of the cloud (and of the mask), projects each valid pixel with `BascOnePoint()`,
then measures the average ground distance to the re-projections of its 8 image
neighbours (`AllocNeighbourhood<2>(1)`) at the same depth/Z. The result is the **median**
over all pixels, truncated to 3 decimals (the raw median is kept when truncation would
yield 0, to avoid a later division by zero). Requires at least 4 valid samples
(asserted). `mBoxLocTarget` is reset afterwards since `BascOnePoint` polluted it.

### `MakeFastBasc()` — mesh construction + Z-buffer (per tile)

1. **`BoxUtile(mImMasq1)`** — tight bounding box of mask pixels `> 0.99`; early return if empty.
2. **Vertex pass** — for every pixel in the useful box, call `BascOnePoint()`. Valid
   pixels get an increasing index in `anImIndex` (`-1` = no vertex, e.g. depth 0) and
   push:
   * `aVPts` — 3D ground point (-> `mTri3D`),
   * `aVPtsDepth` — `(x, y, depth)` point (-> `mTri2DDepth`),
   * `mVCorrel` — correlation value (if `mBascCorrel`).
3. **Face pass** — for each 2×2 pixel block whose 4 corners have valid vertices, create
   two triangles `(I00, I10, I11)` and `(I00, I11, I01)`. Faces are shared *by index*
   between `mTri3D` and `mTri2DDepth` (same `aVFaces`), which is what allows the
   Z-buffer to iterate both meshes in lock-step.
4. **Z-buffer** (see §6) — build `cZBuffer` and run
   `MakeZBufForBasc(ProjInit)` then `MakeZBufForBasc(SurfDevlpt)`.
5. **`ProcessNoPix(aZB)`** — rescue of `NoPix` triangles (see §6.3).
6. **`MakeBasculeTris(aZB)`** — rasterize (see §5.6).

Each stage is wrapped in a `cAutoTimerSegm` so timings appear in the MMVII log
(`ZBuffer::PrepareProjections`, `ZBuffer::Init`, `ZBuffer::ProjInit`,
`ZBuffer::SurfDevlpt`, `ZBuffer::ProcessNoPix`, `ZBuffer::Pixellization`).

Key Z-buffer wiring:

```cpp
// Ground -> (image x, image y, depth-or-Z) mapping, depending on the sensor type
aMapCamDepth = cSensorCamPC ? cSIMap_Ground2ImageAndProf : cSIMap_Ground2ImageAndZ;

cSetVisibility aSetVis(mCamGen, mMII);     // input set  : visibility from the camera
cDataBoundedSet<tREAL8,3> aSetCam(aBox);   // output set : tile box, Z unbounded

cZBuffer aZBuf(aTriIt,        // iterator over mTri3D   (ground mesh)
               aSetVis, *aMapCamDepth, aSetCam,
               /*aResolOut*/ 1,
               /*aSameOri */ true,
               /*IsBascProc*/true,
               aTriIT2DDepth); // iterator over mTri2DDepth (image+depth mesh)
```

Because `IsBascProc=true`, the Z-buffer skips the `Ground->Image` mapping of each
triangle (the *output* triangles are read directly from `mTri2DDepth`) and skips the
resolution computation — the mapping is just depth-map -> DEM.

### `MakeBasculeTris(cZBuffer&)` — rasterization and tile output

* Builds the tile's **ground->pixel affinity**:

  ```cpp
  anAffinetoTarget = cAffin2D(origin = (Xmin, Ymax) of mBoxLocTarget,
                              VX = (GSD, 0), VY = (0, -GSD));   // north-up raster
  ```

  and its inverse `anInvAfftoPixel`. Raster size = `ceil(BoxSz / GSD)`.
* `mImRed` (tREAL8 elevations, initialized to `mInfty = -1e10` = no-data),
  `mImMasqOut`, and optionally `mImCorrelOut` are (re)sized to the tile raster.
* For each face with `ZBufLabIsOk(result)` (i.e. `Visible` or `LikelyVisible`):
  * map the 3 ground vertices to raster pixels, build a `cTriangle2DCompiled`,
  * `PixelsInside(aVPix, 1e-10, &aVW)` gives covered pixels + barycentric weights,
  * interpolated elevation `Scal(aElev, aVW[k])` is written to `mImRed`; the mask is set
    to 1 (last writer wins — occlusions were already resolved by the Z-buffer),
  * with `mBascCorrel`, the vertex correlations are barycentrically interpolated and
    down-weighted by `min(1, 1/mStretchThresh)` of the face; a (near) zero weighted
    correlation zeroes the mask pixel.
* Counters `NB GOOD / NB BAD(Distorted) / NB HIDDEN / NB NO PIX` are printed.
* Tile outputs written into `mPhProj.DPMeshDev().FullDirOut()`, indexed by the tile
  index `mIndBoxRecal`:

  | File | Content | Type |
  |---|---|---|
  | `BLOC-<ix>-<iy>-<NameResult>` | elevation tile | `eTN_REAL8` |
  | `MASQ-<ix>-<iy>-<NameResult>` | validity mask | `eTN_U_INT1` |
  | `CORR-<ix>-<iy>-<NameResult>` | weighted correlation (optional) | `eTN_U_INT1` |
  | `*-<ix>-<iy>-<NameResult .tfw>` | TFW georeferencing of each raster | text |

### `GenTFW` / `ReadTFW` — world-file I/O

A minimalist 6-line TFW writer/reader for a `cAffin2D<tREAL8>`:

```
VX.x  VX.y  VY.x  VY.y  Tr.x  Tr.y     (one value per line, precision 10)
```

`ReadTFW` is `static` (also used by the bench) and asserts the file opens and holds at
least 6 lines. Note the layout is MMVII's affinity ordering, read back as
`cAffin2D(Tr, VX, VY)`.

### `MergeResults()` — global mosaic (master process only)

1. Re-derives the tile grid with `cParseBoxInOut<2>::CreateFromSize(mDFI2d, mSzTiles)`.
2. **Pass 1 — global footprint:** for each tile index that produced a
   `MASQ-...tfw`, read the TFW, compute the useful mask box (`BoxUtile`), map its
   corners to ground and grow `mBoxGlobTarget`. Missing tiles (empty mask) are skipped.
3. Builds the global affinity `mGlobAff` (same north-up convention) and the global
   raster size `mBoxGlobOutPix = ceil(BoxSz / GSD)`.
4. Creates the merged outputs in `mOutDir`:
   * `Prof_<NameResult>` (`eTN_REAL4`, init `-1e9`) + `Prof_....tfw`
   * `Masq_<NameResult>` (`eTN_U_INT1`) + `Masq_....tfw`
   * `Correl_<NameResult>` (+ tfw) when `mBascCorrel`
5. **Pass 2 — mosaicking:** each tile's useful crop is read back and pasted at
   `mGlobAff^-1(localAff(pixel))`. Where tiles overlap, a **max-Z rule** is applied
   (`aCurUpdateZ > aSavedZ` wins), which also discards the `mInfty` no-data values.
6. Rasters are flushed to disk.

Back in `Exe()`, the per-tile `BLOC-* / MASQ-* / CORR-*` files are deleted and a
`cCloudRaster` XML (`<Prof_...>.exml` name pattern `mNameBascOut + ".xml"`) is written
with the merged file names, pixel size and `mGlobAff` — this is the hand-off contract
for `CloudFuser`.

---

## Interaction with the Z-buffer (`cZBuffer`)

Implementation: [`src/Geom3D/cZBuffer.cpp`](../../src/Geom3D/cZBuffer.cpp),
interface: [`include/MMVII_ZBuffer.h`](../../include/MMVII_ZBuffer.h).

### Bascule-specific constructor path (`IsBascProc = true`)

In bascule mode the constructor iterates **all** mesh points without filtering by
`mSetIn`/`mSetOut` (the appli guarantees them valid) and keeps the second mesh
(`mMesh2DDepth`) in sync. The output bounding box, the `Out->Pix` homothety
(`mROut2Pix`, scale `1/mResolOut`, border of 2 px) and the Z-buffer image
(`mZBufIm`, init `mInfty`) are set up as usual.

### `MakeZBufForBasc(eZBufModeIter)`

Called twice by `MakeFastBasc`:

* **`ProjInit`** — first pass; for each (3D, 2D+depth) triangle pair, splat the
  interpolated depth into `mZBufIm`, keeping the max of `mMultZ * depth`
  (`aSameOri=true => mMultZ=-1`, i.e. nearest-to-camera wins).
* **`SurfDevlpt`** — second pass; each triangle is re-rasterized and compared against
  the final Z-buffer:
  * >= 50 % of its pixels at the stored depth -> `Visible`, else `Hidden`;
  * wrong normal orientation -> `BadOriented`;
  * degenerate input/output triangle -> `UnRegIn` / `UnRegOut`;
  * zero covered pixels -> `NoPix`;
  * additionally (bascule only) `IsStretched<2>()` compares the singular values of the
    image->ground 2D deformation; ratio `sigma1/sigma2 > mMaxStretching` (or `sigma2 ~ 0`)
    re-labels the face `Distorted` and records the ratio in
    `cResModeSurfD::mStretchThresh` (used as the correlation down-weight in
    `MakeBasculeTris`).

  One `cResModeSurfD {mResult, mResol, mStretchThresh}` per face is accumulated in
  `mResSurfD`, queried by the appli through `aZB.ResSurfD(kFace)`.

### `ProcessNoPix(cZBuffer&)` — NoPix propagation

Faces smaller than one output pixel get `NoPix` and would be lost. The appli builds the
mesh dual graph (`mTri3D->MakeTopo()`, `cGraphDual`) and iteratively relabels any
`NoPix` face adjacent to an OK face (`ZBufLabIsOk` = `Visible` or `LikelyVisible`) as
**`LikelyVisible`**, looping until a fixed point. This flood-fills visibility across
sub-pixel triangles connected to visible areas while leaving isolated `NoPix` islands
rejected.

---

## Parallelism model

* The master process computes everything global (XML parse, sensors, GSD), then
  `APBI_ExecAll()` spawns one sub-process per 2000×2000 tile (overlap 50 px).
* Values computed once by the master are forwarded to recalls via `mArgsRecall`
  (currently only `GSD=<value>`); each recall therefore sees `IsInit(&mGSD) == true`.
* Each recall executes `Exe()` down to `APBI_ExecAll()` which, inside a recall,
  processes only its own tile (`ExeOnParsedBox`), writes the `BLOC/MASQ/CORR` tile
  files, and returns. Merging, cleanup and the `cCloudRaster` XML are performed only
  when `!InsideParalRecall()`.
* Tiles communicate exclusively through the file system (tile rasters + TFW files in
  `DPMeshDev` out-dir), so there is no shared-memory state to worry about.

---

## Bench

`BenchCloudEpipProjectOnGround(cParamExeBench&)` (same source file, registered in
[`src/Bench/BenchGlob.cpp`](../../src/Bench/BenchGlob.cpp)) validates the EPIP geometry
chain on the synthetic dataset `EPIP/Synthetic/`: it reads the two epipolar cameras and
the disparity map, triangulates a few hand-picked pixel triangles via
`RobustBundleInters`, rasterizes them against the ground-truth DEM georeferencing
(`ReadTFW` of `GroundTruth_Dem_Epip1.tfw`) and checks the interpolated elevations match
the ground-truth raster.
