::: center
**Help for command EpipRectification**
:::

# Description of EpipRectification

This command computes the epipolar geometry of a pair of images and,
optionally, resamples both images (and generates their sensor models) in
this geometry. In epipolar geometry the two images are transformed so
that the rows of the two resampled images match: a ground point has the
same `y` in both images, and only its `x` differs (the disparity). This
is the geometry expected by dense matching.

The command works with any sensor that can project a ground point and
build the bundle of an image point: RPC (satellite) sensors as well as
central perspective cameras. The rectification does not need a regular
epipolar geometry: the transformation of each image is a polynomial
model fitted from synthetic correspondences, generated between the two
sensors over a Z interval.

The result of the computation is a *model* (one polynomial mapping per
image), which can be saved with `SaveModel=true` and reused later by the
command `EpipResampling` (see the last section).

# Principle and conventions

The epipolar mapping of an image is

    (x,y) -> (x_rot, V(x_rot,y_rot)) - frame origin

that is: a rotation about the centre of the image (aligning the epipolar
direction with the `x` axis), followed by a polynomial `V` on `y`
(depending on `x_rot` and `y_rot`), then a translation so that the
coordinates of the output image start at 0. Consequently:

-   `x` is the own rotated coordinate of each image;

-   `y` is common to the pair: the rows match (`y1 = y2`);

-   the disparity `x2 - x1` depends on the ground Z and on the position
    in the image.

The direct polynomials (`V1`, `V2`) are fitted with degree `Degree`; the
inverse polynomials (`W1`, `W2`), needed to resample, with degree
`DegreeInv`. The default of `DegreeInv` is `Degree + 4`.

The correspondences used for the fit are generated from the sensors
themselves: points of the first image are lifted to several Z values of
the Z interval and projected in the second image, and conversely. They
are split into a training pool, used for the fit, and a held-out pool,
used for the quality control (see below). `ZSteps` is the number of Z
steps used for this generation.

The size of the output images is controlled by `FrameAlgo`, which chooses
how the common frame of the two epipolar images is built: `Intersect`
(default) keeps only the common part, `Union` keeps all the parts,
`Img_1` (resp. `Img_2`) takes the height of the frame from the first
(resp. second) image.

# Mandatory parameters

The command has three mandatory parameters, in this order:

-   the name of the first image;

-   the name of the second image;

-   the name of the orientation (sensor) directory of the project. The
    sensor of each image is read with the usual project mechanism: for
    RPC sensors that are not yet in the project, declare them first with
    `ImportInitExtSens`.

# Z interval

The geometry of two images depends on the range of Z of the scene. The
Z interval used is chosen with the following priority:

1.  `ZIntv=[Zmin,Zmax]`, if given: it overrides every other source (a
    warning is issued when it overrides the interval of the sensor);

2.  an interval inferred from tie points, if `TieP=` is given (name of
    the tie point directory);

3.  the interval of the sensor itself (available for an RPC).

When none of these sources exists, for instance for a central
perspective camera without tie points, the command stops with an error:
give `ZIntv=` or `TieP=`.

When the interval is inferred from tie points, the tie points are
triangulated; those whose triangulation residual exceeds `TiePMaxRes`
(pixels) are discarded. The interval is the envelope of the Z of the
retained points, widened by the relative margin `ZMargin` (default
10%). At least `max(TiePMinNbFloor, TiePMinNbRatio * sqrt(W*H))` points
must be retained, otherwise the command stops. For a scene that is
almost flat, the Z interval inferred from tie points can be too narrow:
give `ZIntv=` explicitly.

# Quality control

The fit is controlled on a set of points that were not used to compute
it. The command prints the standard deviation (in pixels) of the
residual of the direct polynomials `V1`, `V2`, and of the inverse
polynomials `W1` and `W2`. If one of them is above `MaxResid` (default
0.1 pixel) the command stops with an error. When it happens, the usual
remedies are to increase `Degree` (and `DegreeInv`), or to check the Z
interval.

# Checks on the pair

Before the fit, the command samples the first image over the Z interval
and counts the points visible in the second one, and conversely. If the
images do not overlap enough (fewer than 50 points in each pool), it
stops with an error that gives the number of points seen and the Z
interval used: check the images, the orientations and `ZIntv`. With less
than 20% of the points visible, a warning reports the small overlap.
After the fit, a warning is issued if an epipolar mapping folds over the
overlap (the sign of its Jacobian changes: the resampled image would be
unusable there); try a lower `Degree` or check the Z interval.

When the info file is written, the command prints the mean parallax over
the Z interval, that is the change of the disparity with Z alone (the
disparity range also contains its variation with the position in the
image). A warning is issued if it is under 1 pixel (nearly parallel
views, or a Z interval too small: no depth information) or over half of
the width of the epipolar image (almost surely a wrong Z interval).

# Resampling and outputs

Unless `NoImage=true` and `NoRPC=true`, both images are resampled with
the interpolator given by `Interpol` (default: `[Cubic,-0.5]`). The
files are written in `OutDir` (default: `VISU/EpipRectification`).

Names are built from the pattern `OutName`, whose default is
`Epip_%1_%2.tif`: `%1` is replaced by the name of the image (without
directory and extension) and `%2` by the name of the other image of the
pair; `.tif` is added if the pattern has no extension (same for
`MaskName`). For the images `Im1.tif` and `Im2.tif`, the defaults give:

-   `Epip_Im1_Im2.tif` and `Epip_Im2_Im1.tif`: the resampled images;

-   `RPC_Epip_Im1_Im2.tif.xml` and `RPC_Epip_Im2_Im1.tif.xml`: the RPC
    sensor of each resampled image, usable in the following steps (not
    generated with `NoRPC=true`);

-   with `SaveModel=true`, the model of the pair, in a single file named
    after the resampled first image, with the suffix `.EpipModel.` and
    the usual tagged suffix of the profile (`Epip_Im1_Im2.EpipModel.xml`
    by default). It also names the RPC files of the two full frames
    (unless `NoRPC=true`), which `EpipResampling` reuses instead of
    fitting them again.

`NoImage=true` suppresses the resampled images (for instance to compute
only the model, or only the RPC), `NoRPC=true` suppresses the RPC.

The validity mask is optional. With `Mask=true`, or when `MaskName=` is
given, a 1-bit image is written for each resampled image: a pixel is 1
when it maps inside the source image, 0 otherwise. `MaskName` is a name
pattern, in which `$1` is replaced by the name of the resampled image
without extension (default `mask_$1.tif`).

# Info file

When the resampled images are written (not with `NoImage=true`), the
command also writes a file named after the resampled image 1, with the
suffix `.Info.` (and the usual tagged suffix). It describes the pair of
images as they are used for dense matching: the names of the two
resampled images, the box `[P0, P1[` of each image in the epipolar
coordinates (here the whole epipolar frame of each image, image 1 being
the master), the size of each resampled image and the size of the
whole epipolar frame of each image, the shift between them, the Z
interval and the *disparity range* (slave minus master) in the pixel
coordinates of the images, i.e.
`d = (x_slave - CropSlave0.x) - (x_master - CropMaster0.x)`. The disparity
range is computed over the whole master frame from the Z interval. It is
the information needed to set the disparity range of a dense matching.
The command prints the Z interval and the disparity range.

`EpipRectification` has no crop: to resample a window of the epipolar
frame, save the model with `SaveModel=true` and use `EpipResampling`,
which has `CropP0`/`CropP1` and also writes this file, with the boxes of
the crops.

The resampled output is held in memory: for very large frames, use the
crop or the tiling of `EpipResampling`.

# Examples

Rectification of a pair of RPC images, saving the model:

    MMVII EpipRectification Im1.tif Im2.tif Ori SaveModel=true

Central perspective cameras, Z interval given:

    MMVII EpipRectification Im1.tif Im2.tif Ori ZIntv=[0,100]

Z interval inferred from the tie points of the project:

    MMVII EpipRectification Im1.tif Im2.tif Ori TieP=Std

With validity masks:

    MMVII EpipRectification Im1.tif Im2.tif Ori SaveModel=true Mask=true

Computation of the model only:

    MMVII EpipRectification Im1.tif Im2.tif Ori \
        SaveModel=true NoImage=true NoRPC=true

# Related commands

`EpipResampling` resamples both images of a pair from the model saved by
`SaveModel=true`, without recomputing the geometry. Only it offers the
crop and the cutting of the pair in tiles resampled in parallel; see its
own help.
