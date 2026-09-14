#ifndef  _MMVII_PLYEXPORT_H_
#define  _MMVII_PLYEXPORT_H_

#include "MMVII_Geom3D.h"     // cPlyVertices
#include "MMVII_AllClassDeclare.h"

namespace MMVII
{
class cSensorImage;           ///< perspective or pushbroom sensor
class cComputeMergeMulTieP;   ///< multi-image tie points, merged
class cPlyVertices;           ///< accumulator of coloured vertices



/* ******************************************** */
/*                                              */
/*            cPlyExportCamGeom                 */
/*                                              */
/* ******************************************** */

/// Export the geometry of a set of cameras (centers, image planes, frustums)
/// as coloured vertices in a ply cloud. Depends only on cSensorImage, so it can
/// be used without any tie point or image data.

class cPlyExportCamGeom
{
public:
    /// aCamScale : multiplies the diagonal/focal ratio, sets pyramid length
    cPlyExportCamGeom(tREAL8 aCamScale,int aNbSteps=50);

    void AddCameras(cPlyVertices& aPlyverts, const std::vector<cSensorImage *>& aVSens) const;

    /// perspectives centers only, in red
    void AddCenters(cPlyVertices & aPlyverts, const std::vector<cSensorImage *>& aVSens) const;
    /// image plane border (+grid if non-pinhole camera), in red
    void AddImagePlanes(cPlyVertices & aPlyverts, const std::vector<cSensorImage *>& aVSens) const;
    /// frustum sampled from center to image plane, in green
    void AddFrustums(cPlyVertices & aPlyverts, const std::vector<cSensorImage *>& aVSens) const;

private:
    /// Depth at which the image plane is drawn, from image size and focal
    tREAL8  CalculateFDepth(const cPt2di & aSz,tREAL8 aF) const;

    tREAL8  mCamScale;   ///< length of the camera pyramid
    int     mNbSteps;    ///< sampling density along image borders and frustum
};

/* ******************************************** */
/*                                              */
/*            cPlyExportTiePoints               */
/*                                              */
/* ******************************************** */

/// Export merged multi-image tie points as a coloured ply cloud : filter the
/// observations by reprojection residual, then read each image once to get RGB.
class cPlyExportTiePoints
{
public:
    cPlyExportTiePoints(tREAL8 aErrProjMax,bool aWithRGB,bool aWithAvgRGB);

    void AddTiePoints
        (
            cPlyVertices &                      aPlyVerts,
            cComputeMergeMulTieP &              aTPts,
            const std::vector<cSensorImage*> &  aVSens
            ) const;

    void SetVerbose(bool aVerbose);   ///< print "(k/n) ImageName" while reading

private:
    tREAL8  mErrProjMax;   ///< max reprojection residual to keep an observation
    bool    mWithRGB;      ///< colorize from images
    bool    mWithAvgRGB;   ///< average RGB over all images instead of one
    bool    mVerbose;

};

};
#endif  //  _MMVII_PLYEXPORT_H_
