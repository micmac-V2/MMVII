#ifndef C_EPIPOLAR_RECTIFICATION_H
#define C_EPIPOLAR_RECTIFICATION_H

#include "cPolyXY_N.h"
#include "MMVII_Mappings.h"
#include "MMVII_Geom3D.h"
#include "MMVII_AllClassDeclare.h"  // cPt2dr, cPt3dr, cPt2di, etc.
#include "MMVII_Tpl_ElemStrToVal.h"
#include "MMVII_MeasuresIm.h"  // cSetHomogCpleIm
#include <optional>
#include <string>

namespace MMVII {

class cSensorImage;
class cSensorCamPC;
class cParamExeBench;

/// Serialization for cRect2 (no existing AddData for cTplBox/cPixBox in MMVII)
void AddData(const cAuxAr2007 &anAux, cRect2 &aRect);

// Common vertical interval for epipolar Image pair resamling
enum class eEpipFrm
{
    eIntersect,         // Contains only common parts
    eUnion,             // Contains all parts
    eImg_1,             // Frame height from Image 1
    eImg_2,             // Frame height from Image 2
    eNbVals
};



class cEpipolarMapping : public cDataInvertibleMapping<tREAL8,2>
{
public:
    cEpipolarMapping(const cPt2dr& aZInterval) : mZInterval(aZInterval) {}
    // Base's copy ctor is deleted but holds no state we rely on ; reconstruct it fresh.
    cEpipolarMapping(const cEpipolarMapping &aOther)
        : mEpipImFrame(aOther.mEpipImFrame), mZInterval(aOther.mZInterval)
    {}
    void SetEpipImFrame(const cRect2& aFrame) { mEpipImFrame = aFrame;}
    cRect2 EpipFrame() const { return mEpipImFrame; }
    cPt2di EpipImSz() const { return mEpipImFrame.Sz(); }
    cPt2dr ZInterval() const { return mZInterval; }

    /// Tag of the type in a serialized pair model
    virtual std::string TypeName() const = 0;
    virtual std::unique_ptr<cEpipolarMapping> Clone() const = 0;
    virtual void AddData(const cAuxAr2007 &anAux) = 0;
    /// Gives the mapping the sensor it applies to when it needs one (after reading a model) ; nothing by default
    virtual void SetSourceSensor(const cSensorImage &) {}
protected:
    /// Serialize the fields common to every cEpipolarMapping ; called by derived classes' own AddData
    void AddDataBase(const cAuxAr2007 &anAux);

    cRect2 mEpipImFrame{cPt2di{0,0},cPt2di{0,0},true}; ///< frame in epipolar space (for resampling)
    cPt2dr mZInterval;
};

// Default mask name pattern : $1 = output image name without extension
static constexpr const char* TheDefaultMaskNamePat = "mask_$1.tif";

// Options shared by EpipRectification and EpipResampling : crop and validity mask
struct cEpipCropMaskOpts
{
    cPt2di      mCropP0{0,0};
    cPt2di      mCropP1{0,0};
    int         mMaster = 1;   ///< Image (1 or 2) in which the crop is given, the other one is derived from Z
    bool        mMask = false;
    std::string mMaskName = TheDefaultMaskNamePat;
    // Set by the command from IsInit (not registered)
    bool        mHasCrop = false;
    bool        mMaskOn = false;
};

class cSensorImage;
class cInterpolator1D;

// Crop of the slave epipolar image matching a master crop over a Z interval : same rows as the master,
// columns from the disparity range, clipped to the slave frame.
struct cEpipSlaveCrop
{
    cPt2di mP0;
    cPt2di mP1;
    cPt2dr mDispRange;   ///< range of x2-x1 (epipolar coordinates, no crop) over the master crop and the Z interval
    tREAL8 mMeanParallax = 0;   ///< mean over the sampled points of the change of x2-x1 over the Z interval (px) : the depth information
};

cEpipSlaveCrop EpipSlaveCrop(const cEpipolarMapping & aMapM, const cEpipolarMapping & aMapS,
                             const cSensorImage & aSIM, const cSensorImage & aSIS,
                             const cPt2di & aMasterP0, const cPt2di & aMasterP1,
                             const cPt2dr & aZIntv, int aMargin=2);

// Description of a pair of crops, written next to the resampled images (for dense matching on the crops).
struct cEpipCropInfo
{
    std::string mNameImMaster;   ///< resampled crop image that carries the given crop (file name)
    std::string mNameImSlave;
    cPt2di      mCropMaster0{0,0}; ///< crops in the epipolar coordinates of each image, [P0,P1[
    cPt2di      mCropMaster1{0,0};
    cPt2di      mCropSlave0{0,0};
    cPt2di      mCropSlave1{0,0};
    cPt2di      mSizeMaster{0,0};  ///< size of the resampled images (crop P1-P0)
    cPt2di      mSizeSlave{0,0};
    cPt2di      mFrameSizeMaster{0,0};  ///< size of the whole epipolar image (frame) of each image
    cPt2di      mFrameSizeSlave{0,0};
    cPt2di      mShift{0,0};     ///< mCropSlave0 - mCropMaster0
    cPt2dr      mZInterval;      ///< Z interval used to derive the slave crop
    /// Disparity range in the pixel coordinates of the cropped images :
    /// d = (x_slave - mCropSlave0.x) - (x_master - mCropMaster0.x)
    cPt2dr      mDispRange;

    void AddData(const cAuxAr2007 &anAux);
    void ToFile(const std::string &aNameFile) const;
    static cEpipCropInfo FromFile(const std::string &aNameFile);
};

void AddData(const cAuxAr2007 &anAux, cEpipCropInfo &anInfo);

cEpipCropInfo MakeEpipCropInfo(const cPt2di & aMasterP0, const cPt2di & aMasterP1,
                               const cEpipSlaveCrop & aSlave, const cPt2dr & aZIntv,
                               const cPt2di & aFrameSizeMaster, const cPt2di & aFrameSizeSlave,
                               const std::string & aNameImMaster, const std::string & aNameImSlave);

// Tiles of size aSz covering [aK0,aK1[ with at least aOverlap between neighbours. All tiles have size aSz,
// the last one being moved back inside the interval; a single tile when the interval is not longer than aSz.
std::vector<std::pair<int,int>> EpipTiles1D(int aK0,int aK1,int aSz,int aOverlap);
/// Same in 2D, over [aP0,aP1[, with an overlap per axis ; tiles in row-major order (x fastest)
std::vector<cRect2> EpipTiles(const cPt2di & aP0,const cPt2di & aP1,const cPt2di & aSz,const cPt2di & aOverlap);

// Name of a tile : "_t<row>_<col>" inserted before the extension of the name pattern (fixed width indices)
std::string EpipTileName(const std::string & aPattern,int aRow,int aCol,int aNbRow,int aNbCol);

// Name pattern of an output image : ".tif" is added when it has no extension (a '.' followed by no '%', '$' or '/')
std::string EpipNameWithExtension(const std::string & aPattern);

// File listing the tiles of an automatic tiling : parameters, and the info file of each tile
struct cEpipTileEntry
{
    int         mRow = 0;
    int         mCol = 0;
    std::string mInfoFile;    ///< cEpipCropInfo file of the tile (name, in the same directory)

    void AddData(const cAuxAr2007 &anAux);
};

struct cEpipTilesInfo
{
    std::string mNameModel;           ///< pair model file
    int         mMaster = 1;          ///< image of the model that carries the tiles (Master= option)
    cPt2di      mSzTiles{0,0};
    cPt2di      mSzOverL{0,0};
    cPt2di      mRegion0{0,0};        ///< tiled region of the master epipolar frame, [P0,P1[
    cPt2di      mRegion1{0,0};
    cPt2dr      mDispRange;           ///< disparity range over all the tiles (cropped coordinates of each tile)
    std::vector<cEpipTileEntry> mTiles;

    void AddData(const cAuxAr2007 &anAux);
    void ToFile(const std::string &aNameFile) const;
    static cEpipTilesInfo FromFile(const std::string &aNameFile);
};

void AddData(const cAuxAr2007 &anAux, cEpipTileEntry &anEntry);
void AddData(const cAuxAr2007 &anAux, cEpipTilesInfo &anInfo);

// Folds of an epipolar mapping on a grid of the master image restricted to its overlap with the slave (mapping extrapolated elsewhere) :
// {number of points where the sign of the Jacobian differs from the majority, number of points tested}
std::pair<int,int> EpipFoldCount(const cEpipolarMapping & aMap, const cSensorImage & aSIM, const cSensorImage & aSIS, const cPt2dr & aZIntv);

// Z interval valid for both images of the pair : intersection of the two, error if empty.
cPt2dr EpipPairZInterval(const cEpipolarMapping & aMap1, const cEpipolarMapping & aMap2);

// Sensor of the full epipolar frame of an image, named aName. Read from aSensorFile if not empty (saved by EpipRectification), else fitted
// from aSI and the Z interval of the mapping. Caller owns it.
cSensorImage * EpipSensor(const cSensorImage & aSI, const cEpipolarMapping & anEpipMap, const std::string & aName, const std::string & aSensorFile = "");

// Resample one image in its epipolar geometry (crop, optional mask, optional sensor). aEpipSI : sensor of the full epipolar frame (EpipSensor), cropped here ; null => none written.
void ResampleEpipImage(const cEpipCropMaskOpts & anOpt,
                       const std::string & aNameIm,
                       const cEpipolarMapping & anEpipMap,
                       const cSensorImage * aEpipSI,
                       const cInterpolator1D & aInterp,
                       bool aNoImage,
                       const std::string & aOutDir,
                       const std::string & aOutBaseName);

// Adds a crop (shift of the origin) on top of a cEpipolarMapping, at resampling time only (never serialized).
class cEpipCropMapping : public cDataInvertibleMapping<tREAL8,2>
{
public:
    cEpipCropMapping(const cEpipolarMapping& aMap, cPt2dr aCropP0)
        : mMap(aMap), mCropP0(aCropP0) {}

    cPt2dr Value(const cPt2dr& aPt) const override { return mMap.Value(aPt) - mCropP0; }
    cPt2dr Inverse(const cPt2dr& aPt) const override { return mMap.Inverse(aPt + mCropP0); }

private:
    const cEpipolarMapping& mMap;
    cPt2dr mCropP0;
};





class cEpipolarModel
{
public:
    virtual ~cEpipolarModel() {}
    const cEpipolarMapping& EpipMap1() const { return (const_cast<cEpipolarModel*>(this))->GetEpipMap1();}
    const cEpipolarMapping& EpipMap2() const { return (const_cast<cEpipolarModel*>(this))->GetEpipMap2();;}
    void ComputeCommonFraming(
        const cTplBox<tREAL8,2> aBox1,
        const cTplBox<tREAL8,2> aBox2,
        eEpipFrm aFrmType=eEpipFrm::eIntersect,
        int aMargin=0
    );
protected:
    virtual cEpipolarMapping& GetEpipMap1() = 0;
    virtual cEpipolarMapping& GetEpipMap2() = 0;

};


template<typename T>
#if __cplusplus >= 202002L
    requires std::derived_from<T, cEpipolarModelBase>
#endif
class cEpipolarModelTpl : public cEpipolarModel
{
public:
    typedef std::unique_ptr<T> Ptr_T;
    cEpipolarModelTpl(Ptr_T aEpipMap1, Ptr_T aEpipMap2)
        : aEpipMap1(std::move(aEpipMap1)),aEpipMap2(std::move(aEpipMap2)  )
    {}

private:
    cEpipolarMapping& GetEpipMap1() override { return *aEpipMap1; }
    cEpipolarMapping& GetEpipMap2() override { return *aEpipMap2; }

    Ptr_T aEpipMap1;
    Ptr_T aEpipMap2;
};



class cEpipPolyMapping: public cEpipolarMapping
{
public:
    cEpipPolyMapping(const cPolyXY_Nd& aV,
                     const cPolyXY_Nd& aW,
                     cPt2dr aCenter,
                     cPt2dr aDir,
                     cPt2dr aZInterval,
                     tREAL8 aGridStep, int aNbStepX, int aNbStepY)
        : cEpipolarMapping(aZInterval)
        , mGridStep(aGridStep), mNbStepX(aNbStepX), mNbStepY(aNbStepY)
        , mV(aV)
        , mW(aW)
        , mCenter{aCenter}
        , mDir{aDir}
    {}
    cEpipPolyMapping(const cEpipPolyMapping &aOther)
        : cEpipolarMapping(aOther)
        , mGridStep(aOther.mGridStep), mNbStepX(aOther.mNbStepX), mNbStepY(aOther.mNbStepY)
        , mV(aOther.mV), mW(aOther.mW), mCenter(aOther.mCenter), mDir(aOther.mDir)
    {}

    cPt2dr Value(const cPt2dr& aPt) const override;
    cPt2dr Inverse(const cPt2dr& aPt) const override;

    /// GenerateData's XY grid : pixel step and step count per axis
    tREAL8 GridStep() const { return mGridStep; }
    int NbStepX() const { return mNbStepX; }
    int NbStepY() const { return mNbStepY; }

    std::string TypeName() const override { return "Poly"; }
    std::unique_ptr<cEpipolarMapping> Clone() const override { return std::make_unique<cEpipPolyMapping>(*this); }
    /// Serialize : base (ZInterval,EpipImFrame) then GridStep,NbStepX/Y,V,W,Center,Dir
    void AddData(const cAuxAr2007 &anAux) override;

private:
    /// (p - C) / D  (complex division = rotation)
    cPt2dr ToRotatedFrame(const cPt2dr& p) const;

    /// q * D + C
    cPt2dr FromRotatedFrame(const cPt2dr& q) const;

    tREAL8 mGridStep;
    int    mNbStepX;
    int    mNbStepY;
    // --- Forward polynomial Vk (image -> epipolar) ---
    cPolyXY_Nd mV;
    // --- Inverse polynomial Wk (epipolar -> image) ---
    cPolyXY_Nd mW;
    cPt2dr mCenter;   ///< centroids of the image point sets
    cPt2dr mDir;      ///< unit epipolar direction per image
};

void AddData(const cAuxAr2007 &anAux, cEpipPolyMapping &aMap);


// Epipolar model of an image pair : both mappings (polynomial or central perspective) + the Ori and image names they were computed from.
class cEpipPairModel
{
public:
    cEpipPairModel() {}   ///< for ReadFromFile
    cEpipPairModel(std::unique_ptr<cEpipolarMapping> aMap1, std::unique_ptr<cEpipolarMapping> aMap2, const std::string &aOriName,
                   const std::string &aNameIm1, const std::string &aNameIm2)
        : mMap1(std::move(aMap1)), mMap2(std::move(aMap2)), mOriName(aOriName), mNameIm1(aNameIm1), mNameIm2(aNameIm2)
    {}
    cEpipPairModel(const cEpipPairModel &aOther);
    cEpipPairModel(cEpipPairModel &&) = default;
    cEpipPairModel & operator = (cEpipPairModel &&) = default;

    /// aNum is 1 or 2
    const cEpipolarMapping &Map(int aNum) const { return (aNum==1) ? *mMap1 : *mMap2; }
    /// Gives the mappings that need it (central perspective) the sensors of the images, to call after reading the model
    void BindSensors(const cSensorImage &aSensor1, const cSensorImage &aSensor2);
    const std::string &ImName(int aNum) const { return (aNum==1) ? mNameIm1 : mNameIm2; }
    const std::string &OriName() const { return mOriName; }
    /// File (name relative to the model file) of the RPC of the full epipolar frame of an image, empty if none was saved
    const std::string &RPCName(int aNum) const { return (aNum==1) ? mRPCName1 : mRPCName2; }
    void SetRPCNames(const std::string &aName1, const std::string &aName2) { mRPCName1 = aName1; mRPCName2 = aName2; }

    void AddData(const cAuxAr2007 &anAux);
    void ToFile(const std::string &aNameFile) const;
    /// Read back a model saved by ToFile()
    static cEpipPairModel FromFile(const std::string &aNameFile);

private:
    std::unique_ptr<cEpipolarMapping> mMap1;
    std::unique_ptr<cEpipolarMapping> mMap2;
    std::string      mOriName;
    std::string      mNameIm1;
    std::string      mNameIm2;
    std::string      mRPCName1;   ///< RPC of the full frame saved with the model (optional)
    std::string      mRPCName2;
};

void AddData(const cAuxAr2007 &anAux, cEpipPairModel &aPairModel);


class cEpipPolyModel : public cEpipolarModelTpl<cEpipPolyMapping>
{
public:
    using cEpipolarModelTpl<cEpipPolyMapping>::cEpipolarModelTpl;
};


// ============================================================
//  Epipolar rectification of a pair of central perspective cameras (cSensorCamPC), closed form.
//  The epipolar image of a camera is the image of a virtual camera with the same centre, a rotation common to the pair
//  (X axis = baseline), a common focal and no distortion : no Z interval, no fit.
//  Epipolar space : pixel of the virtual camera with principal point (0,0), minus the origin of the frame.
// ============================================================

class cEpipMappingPC : public cEpipolarMapping
{
public:
    /// aEpipRot : virtual camera to world ; aZInterval : information only (not used by the mapping)
    cEpipMappingPC(const cSensorCamPC & aCam,const cRotation3D<tREAL8> & aEpipRot,tREAL8 aFocal,const cPt2dr & aZInterval=cPt2dr(0,0));
    /// Placeholder to be filled by AddData, then bound to its camera by SetSourceSensor
    cEpipMappingPC();

    cPt2dr Value(const cPt2dr& aPt) const override;
    cPt2dr Inverse(const cPt2dr& aPt) const override;

    const cRotation3D<tREAL8> & EpipRot() const { return mRot; }
    tREAL8 Focal() const { return mFocal; }
    /// Virtual camera of the full epipolar frame (frame must be set), named aName ; caller owns it, its calibration is deleted at end of application
    cSensorCamPC * VirtualCamera(const std::string & aName) const;

    std::string TypeName() const override { return "PC"; }
    std::unique_ptr<cEpipolarMapping> Clone() const override { return std::make_unique<cEpipMappingPC>(*this); }
    /// Serialize : base then Focal and the three axes of the virtual camera (the source camera is not saved)
    void AddData(const cAuxAr2007 &anAux) override;
    void SetSourceSensor(const cSensorImage & aSensor) override;

private:
    cPt2dr ValueInDomain(const cPt2dr& aPt) const;
    const cSensorCamPC &  Cam() const;
    const cSensorCamPC *  mCam;    ///< Source camera, not owned, null until bound when read from a file
    cRotation3D<tREAL8>   mRot;
    tREAL8                mFocal;
};

class cEpipModelPC : public cEpipolarModelTpl<cEpipMappingPC>
{
public:
    using cEpipolarModelTpl<cEpipMappingPC>::cEpipolarModelTpl;
};



// ============================================================
//  Epipolar rectification for a generic stereo camera pair.
//
//  Reference:
//    Pierrot-Deseilligny & Rupnik,
//    "Epipolar rectification of a generic camera", 2021.
//
//  The method computes two mappings
//
//    Fk : Ik -> Ek  (image -> epipolar space)
//
//  such that  F1(P1).y == F2(P2).y  whenever P1 and P2 are
//  H-compatible (i.e. they could be the projection of the same
//  3D point).
//
//  Each mapping has the form (eq. 23) :
//    Fk(p) = ( p_rot.x ,  Vk(p_rot) )
//    where  p_rot = (p - Ck) / Dk   (complex division = rotation)
//
//  The inverse is (eq. 33) :
//    Gk(e) = Dk * (e.x, Wk(e)) + Ck
// ============================================================

class cEpipolarRectification
{

public:
    // --------------------------------------------------------
    //  User parameters
    // --------------------------------------------------------
    struct cParams
    {
        int      mPolyDegree    = 3;      ///< degree of V polynomials
        int      mPolyDegreeInv = 7;      ///< degree of inverse W polynomials
        int      mNbZLevels     = 3;      ///< number of altitude sampling levels
        eEpipFrm mEpipFrm       = eEpipFrm::eIntersect; ///< Framing type for epipolar images (Resmampling)
        int      mMargin        = 2;      ///< Margin in pixels for epipolar image framing (Resampling)
        std::optional<cPt2dr> mZIntv = std::nullopt; ///< Override Z interval (Zmin,Zmax); mandatory if sensor has none
        std::optional<cSetHomogCpleIm> mHomolPts = std::nullopt; ///< Tie points to infer Zmin/Zmax from (alternative to mZIntv, lower priority)
        tREAL8   mTiePMaxRes     = 2.0;   ///< Max triangulation residual (px) to keep a tie point for Z inference
        tREAL8   mZMargin        = 0.10;  ///< Relative margin added around the Z envelope inferred from mHomolPts
        tREAL8   mTiePMinNbRatio = 0.04;  ///< Min kept-tie-point count = max(mTiePMinNbFloor, mTiePMinNbRatio*sqrt(W*H))
        int      mTiePMinNbFloor = 25;    ///< Absolute floor for the min kept-tie-point count above
        bool     mNoWarnings     = false; ///< Don't generate warnigs: used by Bench
        size_t   mMinNbPairs     = 50;    ///< Min accepted pairs required in each of the train/test pools, per master-camera direction
    };

    // --------------------------------------------------------
    //  Constructor
    // --------------------------------------------------------
    cEpipolarRectification(const cSensorImage& aCam1,
                           const cSensorImage& aCam2,
                           const cParams&      aParams);

    // --------------------------------------------------------
    //  Main entry point
    // --------------------------------------------------------
    cEpipPolyModel Compute();

    int NbPairs12() const { return mNbPairs12; }
    int NbPairs21() const { return mNbPairs21; }
    double V1V2Var() const { return mV1V2Var; }
    double W1Var() const { return mW1Var; }
    double W2Var() const { return mW2Var; }
    /// Independent (held-out) residuals : mean square of V1(q1)-V2(q2), W1(...)-q1.y, W2(...)-q2.y on the test pool
    double V1V2VarIndep() const { return mV1V2VarIndep; }
    double W1VarIndep() const { return mW1VarIndep; }
    double W2VarIndep() const { return mW2VarIndep; }
private:
    // --------------------------------------------------------
    //  Private helper : one H-compatible pair in rotated coords
    // --------------------------------------------------------
    struct cEpiPair
    {
        cPt2dr mP1;   ///< rotated point in I1
        cPt2dr mP2;   ///< rotated point in I2
    };

    // ----------------------------------------------------------
    //  Generate H-compatible pairs (Algorithm 2 of the paper), split into a
    //  train pool (fits V1/V2/W1/W2) and a test pool (EstimateIndepResiduals).
    //  Outputs: aOutPairsTrain/Test pairs, aOutCenterM centroid, aOutDirS direction.
    // ----------------------------------------------------------
    void GenerateData(const cSensorImage &aCamM, const cSensorImage &aCamS,
                      std::vector<cEpiPair> &aOutPairsTrain,
                      std::vector<cEpiPair> &aOutPairsTest,
                      cPt2dr &aOutCenterM,
                      cPt2dr &aOutDirS, cPt2dr &aZInterval,
                      tREAL8 &aOutGridStep, int &aOutNbStepX, int &aOutNbStepY) const;

    // ----------------------------------------------------------
    //  Independent residuals of the fitted V1,V2,W1,W2 on the held-out test pairs.
    // ----------------------------------------------------------
    void EstimateIndepResiduals(
            const std::vector<cEpiPair>& aPairsTest,
            const cPolyXY_Nd& aV1, const cPolyXY_Nd& aV2,
            const cPolyXY_Nd& aW1, const cPolyXY_Nd& aW2);

    /// Memoized tie-point-derived Z interval (see EpipEffectiveZInterval)
    mutable std::optional<cPt2dr> mCachedHomolZIntv;

    // ----------------------------------------------------------
    //  Estimate forward polynomials V1 (with Y-axis identity
    //  constraint) and V2 (unconstrained).
    //
    //  System (eq. 24) :  V1(q1) = V2(q2)
    // ----------------------------------------------------------
    void EstimateForwardPolynomials(
            const std::vector<cEpiPair>& aPairs,
            cPolyXY_Nd&           aV1,
            cPolyXY_Nd&           aV2);

    // ----------------------------------------------------------
    //  Estimate inverse polynomials W1, W2 (eq. 33-34).
    //
    //  System :  Wk( qk.x ,  Vk(qk) ) = qk.y
    // ----------------------------------------------------------
    enum class UseFromPair{PT1,PT2};
    void EstimateInversePolynomial(
            const std::vector<cEpiPair>& aPairs,
            const cPolyXY_Nd&     aVk,
            cPolyXY_Nd&           aWk,
            UseFromPair                  aUsePt);

    // --------------------------------------------------------
    //  Members
    // --------------------------------------------------------
    const cSensorImage& mCam1;
    const cSensorImage& mCam2;
    cParams             mParams;
    int mNbPairs12 = 0; ///< number of H-compatible training pairs from I1 to I2 (for info only)
    int mNbPairs21 = 0; ///< number of H-compatible training pairs from I2 to I1 (for info only)
    double mV1V2Var = 0.0;
    double mW1Var = 0.0;
    double mW2Var = 0.0;
    double mV1V2VarIndep = 0.0;
    double mW1VarIndep = 0.0;
    double mW2VarIndep = 0.0;
};

// Z interval of a master camera aCamM of the pair (aCam1,aCam2) : mZIntv > tie-point-derived > aCamM's own native.
// aCache memoizes the tie-point-derived one (identical for both masters).
cPt2dr EpipEffectiveZInterval(const cEpipolarRectification::cParams & aParams,const cSensorImage & aCamM,
                              const cSensorImage & aCam1,const cSensorImage & aCam2,std::optional<cPt2dr> & aCache);

// Closed-form rectification of a pair of central perspective cameras. Uses the frame, Z interval and tie-point
// parameters of cEpipolarRectification::cParams ; the polynomial ones are ignored.
class cEpipolarRectificationPC
{
public:
    typedef cEpipolarRectification::cParams cParams;
    cEpipolarRectificationPC(const cSensorCamPC & aCam1,const cSensorCamPC & aCam2,const cParams & aParams);
    /// Mappings of both images, common frame set
    cEpipModelPC Compute();
private:
    const cSensorCamPC & mCam1;
    const cSensorCamPC & mCam2;
    cParams              mParams;
};



// Body of the closed form bench, run inside the group BenchEpipolarPC
void BenchEpipolarPCBody(cParamExeBench & aParam);

} // namespace MMVII

#endif // C_EPIPOLAR_RECTIFICATION_H
