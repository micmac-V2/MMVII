::: center
**Help for command EpipResampling**
:::

# Description of EpipResampling

This command resamples the two images of a pair in epipolar geometry,
from the model computed by `EpipRectification` and saved with
`SaveModel=true`. The geometry is not recomputed: the model file holds
the mapping of each image, the names of the two images, the name of
the orientation used to compute it and, for a polynomial model, if they
were computed, the names of the RPC of the full frames. The command has
therefore a single mandatory parameter, the model file.

A model is of one of two types, stored in the file: polynomial (generic
solver, RPC or any sensor) or closed form (two central perspective
cameras, see `EpipRectification`). A closed form model holds only the
parameters of the virtual camera of each image (focal, rotation, frame);
the original cameras are read from the orientation stored in the model.

Splitting the two commands allows to compute the geometry once, then to
produce as many outputs as needed: the full frame, a crop, masks, or a
tiling of the pair (see below). The conventions of the
epipolar geometry (coordinates, disparity, direct and inverse
polynomials) are described in the help of `EpipRectification`.

The images of the model are searched under the names stored in the
model, first in the current directory, then in the project directory.
The sensors are read, with the orientation stored in the model, to
derive the crop of the second image and the disparity range of the info
file (see section Crop).

# Resampling and outputs

Both images are resampled with the interpolator `Interpol` (default
`[Cubic,-0.5]`). Files are written in `OutDir` (default
`VISU/EpipResampling`), with names built from the pattern `OutName`
(default `Epip_%1_%2.tif`, where `%1` is the image and `%2` the other
image of the pair, without directory nor extension; `.tif` is added if
the pattern, or `MaskName`, has no extension): the resampled
images, and for each of them the RPC of the resampled image
(`RPC_Epip_Im1_Im2.tif.xml` with the default pattern, unless
`NoOri=true`). The RPC of the full frames saved by `EpipRectification`
are not computed again: their names are stored in the model, the files
are searched next to it, and they are cropped if there is a crop. If one
is missing (model saved with `NoOri=true`, file removed), the RPC is
fitted again, with a warning.

With a closed form model there is no RPC: the sensor of a resampled image
is its virtual camera, rebuilt from the model (no fit, nothing to read)
and written as a standard orientation, `Ori-PerspCentral-Epip_Im1_Im2.tif.xml`
with the default pattern, and its calibration file, in `OutDir`.
`NoOri=true` suppresses it. A crop shifts the principal point of the
camera, the other parameters being unchanged.

With `Mask=true`, or when `MaskName=` is given, a 1-bit mask of the
valid pixels is also written for each image (a pixel is valid when it
maps inside the source image); `$1` in `MaskName` stands for the name of
the resampled image without extension (default `mask_$1.tif`).

# Crop

`CropP0` and `CropP1`, given together, restrict the output to the
half-open box `[CropP0, CropP1[` of the epipolar coordinates of the
image selected by `Master` (1 or 2, default 1). The crop of the other
image is derived from the Z interval of the model. Without a crop, the
whole frames are resampled. In both cases the command writes the info
file of the pair (suffix `.Info.`, named after the resampled image
selected by `Master`), with the boxes of the crops and the disparity
range, as `EpipRectification` does (see its help; without a crop the
boxes are the whole frames). The RPC of a crop shares the polynomials of
the RPC of the full frame, and differs only by its image offsets; the
camera of a crop of a closed form model differs from the full frame one
only by its principal point.

Only the window of the source image needed by a crop is read from the
file.

# Tiling

For images too big to be resampled in one piece, or to feed a dense
matching tile by tile, `SzTiles=[Sx,Sy]` cuts the pair in tiles, which
are resampled in parallel (the number of processes is `NbProc`).

-   the region that is tiled is the whole epipolar frame of the image
    selected by `Master`, or the box given by `CropP0`/`CropP1`;

-   all tiles have the size `SzTiles`; neighbouring tiles overlap by at
    least `SzOverL` (default 0, must be smaller than `SzTiles`); the
    last tile of a row or column is moved back inside the region, so it
    overlaps more rather than being smaller;

-   each tile is a crop of the pair: the command calls itself once per
    tile, with the crop of the tile, and the other arguments unchanged.
    A tile never cuts again;

-   the tiles are named from `OutName`, by inserting `_tRR_CC` (row and
    column, with fixed width) before the extension, so that with the
    default pattern the tile of row 1, column 2 of a 3 by 3 tiling of
    `Im1` and `Im2` is `Epip_Im1_Im2_t01_02.tif`.

Each tile writes its own file of crop and disparity range. When all the
tiles are done, the command writes an *index* of the tiling, in a file
named after the resampled master image, with the suffix `.EpipTiles.`
(and the usual tagged suffix). It gives the model, the master image, the
size and overlap of the tiles, the region, the file of each tile (row,
column), and the disparity range of the whole tiling, which is the union
of those of the tiles. If one tile fails, or has no result, no index is
written and the command stops with an error.

# Examples

Resampling of the pair from the model saved by `EpipRectification`:

    MMVII EpipResampling \
        VISU/EpipRectification/Epip_Im1_Im2.EpipModel.xml

Crop with masks:

    MMVII EpipResampling \
        VISU/EpipRectification/Epip_Im1_Im2.EpipModel.xml \
        CropP0=[0,0] CropP1=[2000,1500] Mask=true

Tiling, with 4 processes:

    MMVII EpipResampling \
        VISU/EpipRectification/Epip_Im1_Im2.EpipModel.xml \
        SzTiles=[2000,2000] SzOverL=[100,100] NbProc=4

# Related commands

`EpipRectification` computes the epipolar geometry of a pair and saves
the model used by this command. It can also resample the images itself,
in one step, with the same mask options but no crop.
