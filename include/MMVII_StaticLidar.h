#ifndef  _MMVII_STATICLIDAR_H_
#define  _MMVII_STATICLIDAR_H_

#include <functional>

#include "MMVII_2Include_Serial_Tpl.h"
#include "MMVII_Geom3D.h"
#include "MMVII_PCSens.h"
#include <unordered_set>

namespace MMVII
{

/** \file MMVII_StaticLidar.h
    \brief Static lidar internal representation

*/

tREAL8 toMinusPiPlusPi(tREAL8 aAng, tREAL8 aOffset = 0.);


class cStaticLidarImporter
{
    friend class cAppli_ImportTSL;
public:
    cStaticLidarImporter();
    void readPlyPoints(std::string aPlyFileName, bool aForceGreenAsIntensity);
    void readE57Points(std::string aE57FileName, bool aForceGreenAsIntensity);
    void readPtxPoints(std::string aPtxFileName, bool aForceGreenAsIntensity);

    /// Adapts to adequate function from postfix, return if some read suceeded
    bool read(const std::string & aName, bool OkNone=false, bool aForceStructured=false,
              std::string aStrInput2TSL="ijk", bool aForceGreenAsIntensity=false);

    void convertToThetaPhiDist();
    void convertToXYZ();
    void ComputeRotInput2Raster(std::string aTransfoIJK); //< give frame used to return to primary rotation axis = z

    bool HasCartesian() const {return mHasCartesian;}
    bool HasIntensity() const {return mHasIntensity;}
    bool HasSpherical() const {return mHasSpherical;}
    bool HasRowCol() const {return mHasRowCol;}
    bool NoMiss() const {return mNoMiss;}
    bool AllPointsReturn() const {return mAllPointsReturn;}
    bool IsStructured() const {return mIsStrucured;}
    int NbCol() const {return mNbCol;}
    int NbLine() const {return mNbLine;}
    tREAL8 ThetaStart() const {return mThetaStart;}
    tREAL8 ThetaStep() const {return mThetaStep;}
    tREAL8 PhiStart() const {return mPhiStart;}
    tREAL8 PhiStep() const {return mPhiStep;}
    tREAL8 DistMinToExist() const {return mDistMinToExist;}
    const std::optional<tPoseR> & ReadPose() const { return mReadPose;}
    bool checkLineCol(); // verify that mMaxCol/mMaxLine ar compatible with mVectPtsLine/mVectPtsCol
    void decimXY(const cPt2di & aDecim);
    const cRotation3D<tREAL8> & RotInput2TSL() const { return mRotInput2TSL; }
    const cRotation3D<tREAL8> & RotInput2Raster() const { return mRotInput2Raster; }

    float ColToLocalThetaApprox(float aCol) const;
    float LineToLocalPhiApprox(float aLine) const;
    float LocalThetaToColApprox(float aTheta) const;
    float LocalPhiToLineApprox(float aPhi) const;
    void ComputeAgregatedAngles();
    float LocalPhiToLinePrecise(float aPhi) const;
    float LocalThetaToColPrecise(float aTheta) const;
    cPt2dr Input3DtoRasterAngle(const cPt3dr & aPt3DInput) const;
    static std::string DefaultPoseName(const std::string & aDirStaticLidarRasters, const std::string & aLidarId);

    void MakeIdImage(const std::string & aNameFile) const; // create miniature to select this scan along images

    // line and col for each point
    std::vector<int> mVectPtsLine;
    std::vector<int> mVectPtsCol;
    // points
    std::vector<cPt3dr> mVectPtsXYZ;
    std::vector<tREAL8> mVectPtsIntens; // 0-1, may be empty if no intensity
    std::vector<cPt3dr> mVectPtsTPD;

    // agregated angles per col/line
    std::vector<tREAL8> mVectPhisCol;
    std::vector<tREAL8> mVectThetasLine;

protected:
    // data
    bool mHasCartesian; //< in original read data
    bool mHasIntensity; //< in original read data
    bool mHasSpherical; //< in original read data
    bool mHasRowCol;    //< in original read data

    bool mNoMiss; // seems to be full (even if some points are (0,0,0)
    bool mAllPointsReturn; // some points are (0,0,0) => no angle!
    bool mIsStrucured;
    std::optional<tPoseR> mReadPose;
    tREAL8 mDistMinToExist;

    int mNbCol, mNbLine;
    tREAL8 mThetaStart, mThetaStep;
    tREAL8 mPhiStart, mPhiStep;
    cRotation3D<tREAL8> mRotInput2TSL; // from xyz input file to classical TSL (rot around z)
    cRotation3D<tREAL8> mVertRot; //< verticalizarion rotation in cloud frame
    cRotation3D<tREAL8> mRotInput2Raster; //< to go from z vertical to z view direction of PP, and make PPx in center
};

// record all data for each patch
struct cLidarRasterPatch
{
    size_t                          mId;        //< Number in cStaticLidar.mPatchCenters
    std::vector<cPt2di>             mLPatchesP; //< px in raster, consituted by points in a lidar scan, begin() is center
    std::unordered_set<std::string> mHiddenOnImage; //< for Im/scanB names, if hidden via zbuffers
    cPt3dr                          mNormalInstr; //< normal in instrument frame
};


class cStaticLidar: public cSensorCamPC
{
    friend class cAppli_ImportTSL;
    friend class cStaticLidarImporter;
public :

    cStaticLidar(const std::string &aNameFile, const std::string & aStationName,
                 const std::string & aScanName, const tPose &aPose, cPerspCamIntrCalib *aCalib,
                 cRotation3D<tREAL8> aRotInput2Raster, tREAL8 aSigma);
    ~cStaticLidar();

    static cStaticLidar *FromFile(const std::string & aNameScanFile, bool aSVP);
    bool ReadRasters(const std::string & aDataDir);

    bool AreRastersReady() const { return mAreRastersReady;}
    void ToPly(const std::string & aName, bool useMask=false) const;
    void AddData(const  cAuxAr2007 & anAux) ;
    virtual void ToFile(const std::string &) const override;
    static std::string RasterIntensityPath(const std::string & aImName); ///< base name should be Station-scan
    static std::string RasterIntensityPath(const cPhotogrammetricProject & aPhProj, const std::string & aImIDName); ///< base name should be Station-scan
    void FillRasters(const cStaticLidarImporter & aSL_importer);
    void SaveRasters(const cStaticLidarImporter & aSL_importer, const std::string &aPhProjDirOut);
    static std::string NameFromId(const std::string &aIdName, bool getOriName);
    static bool IsNameTSL(const std::string &aImageName);

    cCalculator<double> * CreateEqColinearity(bool WithDerives, int aSzBuf, bool ReUse) override; // colinearity with fixed F and PP
    cCalculator<double> * CreateEqColinearityDist(bool WithDerives, int aSzBuf, bool ReUse);
    cCalculator<double> * GetEqColinearityDist();
    void PushOwnObsColinearity(std::vector<double> & aVObs, const cPt3dr &) override; // use this for GCP obs
    void PushOwnObsColinearityDistance(std::vector<double> & aVObs, tREAL4 aMesDistance); // use this for GCP obs


    //inline tREAL8 lToPhiApprox(int l, double aPhiStart, double aPhiStep) const { return aPhiStart + l * aPhiStep; }
    //inline tREAL8 cToThetaApprox(int c, double aThetaStart, double aThetaStep) const { return aThetaStart + c * aThetaStep; }

    void FilterIntensity(const cStaticLidarImporter & aSL_importer, tREAL8 aLowest, tREAL8 aHighest); // add to mRasterMask
    void FilterIncidence(const cStaticLidarImporter &aSL_importer, tREAL8 aAngMax);
    void FilterDistance(tREAL8 aDistMin, tREAL8 aDistMax);
    void MaskBuffer(const cStaticLidarImporter &aSL_importer, tREAL8 aAngBuffer, const std::string &aPhProjDirOut);
    void SelectPatchCenters1(int aNbPatches);
    void SelectPatchCenters2(int aNbPatches, cDataIm2D<tU_INT1> *aSupMaskDIm=nullptr);
    void SelectPatchCenters3(int aNbPatches, cDataIm2D<tU_INT1> *aSupMaskDIm=nullptr);
    void MakeVisu(const cPhotogrammetricProject & aPhProj) const;     ///< show 8bit dist image with patch centers
    void MakePatches(std::list<cLidarRasterPatch> &aLPatches, const std::vector<cSensorCamPC *> &aVCam, int aNbPointByPatch, int aSzMin) const;
    std::tuple<tREAL8,tREAL8,tREAL8> AvgDistNbValidAndNbNotMasked() const; //< return average dist for valid points, number of valid points and number of not-masked points

    cPt3dr Image2InputXYZ(cPt2di aRasterPxI) const; // in input frame
    cPt3dr Image2InputXYZ(cPt2dr aRasterPx) const; // in input frame
    cPt3dr Image2InputXYZ_InterpoleDist(cPt2dr aRasterPx, const cInterpolator1D *anInterpol) const; // just fix distance with interpolator on dist image

    template <typename TYPE>
    cPt3dr Image2Camera3D(const TYPE & aRasterPx) const; // in sensor frame (Z forward)
    cPt3dr Image2Camera3D_InterpoleDist(const cPt2dr & aRasterPx, const cInterpolator1D *anInterpol) const; // in sensor frame (Z forward)
    cPt3dr Image2NormalInstr(const cPt2dr &aRasterPx) const; //< normal in sensor frame

    template <typename TYPE>
        cPt3dr Image2ThetaPhiDist(const TYPE & aRasterPx) const;

    cPt3dr Image2Ground(const cPt2di &aRasterPxI) const;
    cPt3dr Image2Ground(cPt2dr aRasterPx) const;
    cPt3dr Image2Ground_InterpoleDist(cPt2dr aRasterPx, const cInterpolator1D *anInterpol) const;

    tREAL4 Image2Distance(cPt2dr aRasterPx) const;
    tREAL4 Image2DistanceInterpol(cPt2dr aRasterPx, const cInterpolator1D *anInterpol) const;
    cPt3dr ImageAndDepth2Ground(const cPt3dr & ) const override;

    cPt2dr Ground2Image(const cPt3dr &aGroundPt) const override;
    cPt3dr Ground2ImageAndDepth(const cPt3dr &) const override;

    void FixPtPxLoopAroundPP(cPt2dr &aPtPx) const override;

    void TriangulateRegular(const std::string &aVisuPath, int aFactor=16);
    void Triangulate(const std::string &aVisuPath, int aFactor=16);
    cTriangulation3D<tREAL8> * getTriangulation() const;

    static std::string  PrefixName() ;
    std::string  V_PrefixName() const override;
    static std::string Pat2Sup(const std::string & aPatSelect);

    cDataIm2D<tREAL4> &getRasterDistance() const;
    bool IsValidPoint(const cPt2dr &aRasterPx) const; ///< is dist>0
    bool IsValidPoint(const cPt2di &aRasterPx) const; ///< is dist>0
    bool IsMaskedPoint(const cPt2dr &aRasterPx) const;
    tREAL8 Sigma() const;
    const std::vector<cPt2di> & PatchCenters() const;

    bool isInsideNormalInterpolator(cPt2dr & aPt) const;
    cDiffInterpolator1D * getNormalInterpolator(const cPt2dr &aPt) const;

    void Show() const override;
    void initInterpolators();
    static std::string GetIdSuffix();
    static std::string GetIdSuffixRegex();

    cIm2D<tU_INT1> projectIntensityFrom(const cStaticLidar& aFrom) const;

    virtual bool DoAddCalibToUk() const override;

    std::tuple<double, double, cPt3dr> getDistSigmaNormalPlane(cPt2dr aCenter, const cPixBox<2> &aPixBox) const; ///< Adjust a plane on defined points

private :
    template <typename TYPE> static void fillRaster(const cStaticLidarImporter & aSL_importer,
                    std::function<TYPE (int)> func, std::unique_ptr<cIm2D<TYPE>> & aIm); // keep image in memory

    cPt2dr Ground2ImagePrecise(const cPt3dr & aGroundPt) const;

    std::string mStationName;
    std::string mScanName;
    std::string mRasterDistancePath;
    std::unique_ptr<cIm2D<tREAL4>> mRasterDistance;
    std::string mRasterIntensityPath;
    std::unique_ptr<cIm2D<tU_INT1>> mRasterIntensity;
    std::string mRasterMaskPath;
    std::unique_ptr<cIm2D<tU_INT1>> mRasterMask;
    std::string mRasterXPath;
    std::unique_ptr<cIm2D<tREAL4>> mRasterX;
    std::string mRasterYPath;
    std::unique_ptr<cIm2D<tREAL4>> mRasterY;
    std::string mRasterZPath;
    std::unique_ptr<cIm2D<tREAL4>> mRasterZ;
    //std::string mRasterThetaPath;
    //std::string mRasterPhiPath;
    //std::string mRasterThetaErrPath;
    //std::string mRasterPhiErrPath;

    bool mAreRastersReady;

    tREAL8 mSigma;   ///< a priori precision on instrument distances
    std::vector<cPt2di> mPatchCenters;

    // rasters for filtering
    std::unique_ptr<cIm2D<tU_INT1>> mRasterMaskBuffer;
    std::unique_ptr<cIm2D<tREAL4>> mRasterScore; // updated on each filter, used to find patch centers. High=bad

    cRotation3D<tREAL8> mRotInput2Raster; //< to go from z vertical to z view direction of PP, and make PPx in center
    // triangulation for patches selection
    cTriangulation3D<tREAL8> * mTriangulation; ///< triangulation of the raster, for zbuffer

    std::vector<std::pair<double,cDiffInterpolator1D *>> mVInterpN;         ///< Interpolators, used to get normal of the images, with min dist to use them (must be sorted by dist)
    cCalculator<double> * mEqDistColinearityDist;
};

template <typename TYPE>
cPt3dr cStaticLidar::Image2Camera3D(const TYPE & aRasterPx) const
{
    cPt3dr aPtInput3D = Image2InputXYZ(aRasterPx);
    cPt3dr aPtCam3D = mRotInput2Raster.Value(aPtInput3D);
    //std::cout<<"   Image > TLS > Camera3D: " << aRasterPx << " => "<< aPtTLS3D <<" => "<< aPtCam3D <<"\n";
    return aPtCam3D;
}

template <typename TYPE>
    cPt3dr cStaticLidar::Image2ThetaPhiDist(const TYPE & aRasterPx) const
{
    cPt3dr aPtCam3D = Image2Camera3D(aRasterPx);
    tREAL8 aDist = Norm2(aPtCam3D);
    if (aDist==0)
        return {0.,0.,0.}; // InternalCalib()->Value() will make an error
    cPt2dr aPx = InternalCalib()->Value(aPtCam3D);
    cPt2dr aDir = (aPx - InternalCalib()->PP()) / InternalCalib()->F();
    if (aDir.x()<-M_PI)
        aDir.x()+=2*M_PI;
    if (aDir.x()>M_PI)
        aDir.x()-=2*M_PI;
    //std::cout<<"Image2ThetaPhiDist "<<aRasterPx<<" "<<aPtCam3D
    //          <<" "<<aPx<<" "<<InternalCalib()->PP()<<" "<<aDir<<"\n";
    return {aDir.x(), aDir.y(), aDist};
}

void AddData(const  cAuxAr2007 & anAux,cStaticLidar & aSL);

}

#endif  //  _MMVII_STATICLIDAR_H_
