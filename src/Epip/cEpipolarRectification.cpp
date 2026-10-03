#include "cEpipolarRectification.h"
#include "MMVII_Error.h"
#include "MMVII_Geom2D.h"
#include "MMVII_Sensor.h"
#include "MMVII_PCSens.h"
#include "MMVII_Interpolators.h"
#include "../Sensors/cExternalSensor.h"

#include <cmath>
#include <cassert>
#include <algorithm>

namespace MMVII {

// Guard of the RPC fit (cRPC.cpp), read and changed by the benches
extern std::vector<std::pair<int,int>> TheRPCFitGrids;
extern tREAL8 TheRPCFitMaxResPx;
extern int TheRPCFitNbRetry;
extern tREAL8 TheRPCFitLastMaxRes;

template<> cE2Str<eEpipFrm>::tMapE2Str cE2Str<eEpipFrm>::mE2S
    {
     {eEpipFrm::eIntersect,"Intersect"},
     {eEpipFrm::eUnion,"Union"},
     {eEpipFrm::eImg_1,"Img_1"},
     {eEpipFrm::eImg_2,"Img_2"},
     };

MACRO_INSTANTIATE_STRIO_ENUM(eEpipFrm,"EpipFrame")


cPt2dr cEpipPolyMapping::ToRotatedFrame(const cPt2dr &p) const
{
    return (p - mCenter) / mDir;
}

cPt2dr cEpipPolyMapping::FromRotatedFrame(const cPt2dr& q) const
{
    return q * mDir + mCenter;
}



cPt2dr cEpipPolyMapping::Value(const cPt2dr& aPt) const
{
    const cPt2dr q = ToRotatedFrame(aPt);
    return cPt2dr(q.x(), mV.Eval(q)) - ToR(mEpipImFrame.P0());
}

cPt2dr cEpipPolyMapping::Inverse(const cPt2dr& aPt) const
{
    auto aPtEpip = aPt + ToR(mEpipImFrame.P0());
    return FromRotatedFrame(cPt2dr(aPtEpip.x(), mW.Eval(aPtEpip)));
}

// ============================================================
//  Serialization
// ============================================================

void AddData(const cAuxAr2007 &anAux, cRect2 &aRect)
{
    cPt2di aP0 = aRect.P0();
    cPt2di aP1 = aRect.P1();
    MMVII::AddData(cAuxAr2007("P0", anAux), aP0);
    MMVII::AddData(cAuxAr2007("P1", anAux), aP1);
    if (anAux.Input())
        aRect = cRect2(aP0, aP1, true);
}

void cEpipolarMapping::AddDataBase(const cAuxAr2007 &anAux)
{
    MMVII::AddData(cAuxAr2007("ZInterval", anAux), mZInterval);
    MMVII::AddData(cAuxAr2007("EpipImFrame", anAux), mEpipImFrame);
}

void cEpipPolyMapping::AddData(const cAuxAr2007 &anAux)
{
    AddDataBase(anAux);
    MMVII::AddData(cAuxAr2007("GridStep", anAux), mGridStep);
    MMVII::AddData(cAuxAr2007("NbStepX", anAux), mNbStepX);
    MMVII::AddData(cAuxAr2007("NbStepY", anAux), mNbStepY);
    MMVII::AddData(cAuxAr2007("V", anAux), mV);
    MMVII::AddData(cAuxAr2007("W", anAux), mW);
    MMVII::AddData(cAuxAr2007("Center", anAux), mCenter);
    MMVII::AddData(cAuxAr2007("Dir", anAux), mDir);
}

void AddData(const cAuxAr2007 &anAux, cEpipPolyMapping &aMap)
{
    aMap.AddData(anAux);
}

// A mapping is saved with a tag of its type, read back through this factory
static void AddDataMapping(const cAuxAr2007 &anAux, const std::string &aTag, std::unique_ptr<cEpipolarMapping> &aMap)
{
    cAuxAr2007 aAuxMap(aTag, anAux);
    std::string aType = anAux.Input() ? std::string() : aMap->TypeName();
    MMVII::AddData(cAuxAr2007("Type", aAuxMap), aType);
    if (anAux.Input())
    {
        if (aType == "Poly")
            aMap = std::make_unique<cEpipPolyMapping>(cPolyXY_Nd(1), cPolyXY_Nd(1), cPt2dr(0,0), cPt2dr(1,0), cPt2dr(0,0), 1.0, 1, 1);   // placeholder, AddData overwrites it
        else if (aType == "PC")
            aMap = std::make_unique<cEpipMappingPC>();
        else
            MMVII_UserError(eTyUEr::eOpenFile, "Unknown type of epipolar mapping : " + aType);
    }
    aMap->AddData(aAuxMap);
}

void cEpipPairModel::AddData(const cAuxAr2007 &anAux)
{
    AddDataMapping(anAux, "Mapping1", mMap1);
    AddDataMapping(anAux, "Mapping2", mMap2);
    MMVII::AddData(cAuxAr2007("OriName", anAux), mOriName);
    MMVII::AddData(cAuxAr2007("Image1", anAux), mNameIm1);
    MMVII::AddData(cAuxAr2007("Image2", anAux), mNameIm2);
    MMVII::AddData(anAux, "RPCName1", mRPCName1, std::string());   // absent when no RPC was saved
    MMVII::AddData(anAux, "RPCName2", mRPCName2, std::string());
}

cEpipPairModel::cEpipPairModel(const cEpipPairModel &aOther)
    : mMap1(aOther.mMap1->Clone()), mMap2(aOther.mMap2->Clone())
    , mOriName(aOther.mOriName), mNameIm1(aOther.mNameIm1), mNameIm2(aOther.mNameIm2)
    , mRPCName1(aOther.mRPCName1), mRPCName2(aOther.mRPCName2)
{
}

void cEpipPairModel::BindSensors(const cSensorImage &aSensor1, const cSensorImage &aSensor2)
{
    mMap1->SetSourceSensor(aSensor1);
    mMap2->SetSourceSensor(aSensor2);
}

void AddData(const cAuxAr2007 &anAux, cEpipPairModel &aPairModel)
{
    aPairModel.AddData(anAux);
}

void cEpipPairModel::ToFile(const std::string &aNameFile) const
{
    // Non-fixed formatting avoids std::fixed truncating tiny coefficients to 0.
    PushPrecTxtSerial(std::numeric_limits<tREAL8>::max_digits10);
    SetFixedFloatTxtSerial(false);
    SaveInFile(*this, aNameFile);
    SetFixedFloatTxtSerial(true);
    PopPrecTxtSerial();
}

cEpipPairModel cEpipPairModel::FromFile(const std::string &aNameFile)
{
    cEpipPairModel aPairModel;
    ReadFromFile(aPairModel, aNameFile);
    return aPairModel;
}


// ============================================================
//  Slave crop of a master crop
// ============================================================

cPt2dr EpipPairZInterval(const cEpipolarMapping & aMap1, const cEpipolarMapping & aMap2)
{
    const cPt2dr aZ1 = aMap1.ZInterval();
    const cPt2dr aZ2 = aMap2.ZInterval();
    const cPt2dr aRes(std::max(aZ1.x(),aZ2.x()),std::min(aZ1.y(),aZ2.y()));
    MMVII_INTERNAL_ASSERT_User(aRes.x() < aRes.y(), eTyUEr::eUnClassedError,
        "The Z intervals of the two images of the pair do not intersect : " + ToStr(aZ1) + " and " + ToStr(aZ2));
    return aRes;
}

cEpipSlaveCrop EpipSlaveCrop(const cEpipolarMapping & aMapM, const cEpipolarMapping & aMapS,
                             const cSensorImage & aSIM, const cSensorImage & aSIS,
                             const cPt2di & aMasterP0, const cPt2di & aMasterP1,
                             const cPt2dr & aZIntv, int aMargin)
{
    // Master points of the crop (corners included) -> ground at several Z -> slave epipolar
    const int aNbSamp = 9;
    cPt2dr aXRange(1e30,-1e30), aDRange(1e30,-1e30);
    int aNbOk = 0;
    tREAL8 aSumParallax = 0;
    int aNbParallax = 0;
    for (int aKx=0 ; aKx<aNbSamp ; aKx++)
        for (int aKy=0 ; aKy<aNbSamp ; aKy++)
        {
            // upper bounds are excluded pixels : last sample on the last pixel
            const cPt2dr aQM(aMasterP0.x() + (aMasterP1.x()-1-aMasterP0.x()) * aKx / (aNbSamp-1.0),
                             aMasterP0.y() + (aMasterP1.y()-1-aMasterP0.y()) * aKy / (aNbSamp-1.0));
            const cPt2dr aPM = aMapM.Inverse(aQM);
            if (aSIM.PixelDomain().Insideness(aPM) <= 0)   // outside the master image : masked, no counterpart needed
                continue;
            tREAL8 aDMin = 1e30, aDMax = -1e30;   // disparity range of this point over Z
            for (const tREAL8 aZ : {aZIntv.x(), (aZIntv.x()+aZIntv.y())/2.0, aZIntv.y()})
            {
                const cPt2dr aQS = aMapS.Value(aSIS.Ground2Image(aSIM.ImageAndZ2Ground(TP3z(aPM,aZ))));
                aXRange = cPt2dr(std::min(aXRange.x(),aQS.x()),std::max(aXRange.y(),aQS.x()));
                const tREAL8 aDisp = aQS.x() - aQM.x();
                aDRange = cPt2dr(std::min(aDRange.x(),aDisp),std::max(aDRange.y(),aDisp));
                aDMin = std::min(aDMin,aDisp);
                aDMax = std::max(aDMax,aDisp);
                aNbOk++;
            }
            aSumParallax += aDMax - aDMin;
            aNbParallax++;
        }
    MMVII_INTERNAL_ASSERT_User(aNbOk>0, eTyUEr::eUnClassedError, "The master region " + ToStr(aMasterP0) + " " + ToStr(aMasterP1) + " has no point inside the master image");

    const cPt2di aSzS = aMapS.EpipImSz();
    cEpipSlaveCrop aRes;
    aRes.mP0 = cPt2di(std::max(0,round_down(aXRange.x())-aMargin), aMasterP0.y());
    aRes.mP1 = cPt2di(std::min(aSzS.x(),round_up(aXRange.y())+1+aMargin), aMasterP1.y());
    aRes.mDispRange = aDRange;
    aRes.mMeanParallax = aSumParallax / std::max(1,aNbParallax);
    MMVII_INTERNAL_ASSERT_User(aRes.mP0.x() < aRes.mP1.x(), eTyUEr::eUnClassedError,
        "The master region " + ToStr(aMasterP0) + " " + ToStr(aMasterP1) + " has no counterpart in the other image for the Z interval "
        + ToStr(aZIntv) + " (do the two images overlap over this Z interval ?)");
    return aRes;
}


// ============================================================
//  Tiles of a region
// ============================================================

std::vector<std::pair<int,int>> EpipTiles1D(int aK0,int aK1,int aSz,int aOverlap)
{
    MMVII_INTERNAL_ASSERT_User(aK1>aK0, eTyUEr::eUnClassedError, "Empty region to tile");
    MMVII_INTERNAL_ASSERT_User((aOverlap>=0) && (aSz>aOverlap), eTyUEr::eUnClassedError,
        "Tile size must be larger than the overlap (non negative) : size " + ToStr(aSz) + ", overlap " + ToStr(aOverlap));
    if (aK1-aK0 <= aSz)
        return {{aK0,aK1}};

    std::vector<std::pair<int,int>> aRes;
    const int aStep = aSz - aOverlap;
    for (int aK=aK0 ; ; aK+=aStep)
    {
        if (aK+aSz >= aK1)     // last tile : aligned on the end, so its overlap with the previous one is >= aOverlap
        {
            aRes.push_back({aK1-aSz,aK1});
            break;
        }
        aRes.push_back({aK,aK+aSz});
    }
    return aRes;
}

std::vector<cRect2> EpipTiles(const cPt2di & aP0,const cPt2di & aP1,const cPt2di & aSz,const cPt2di & aOverlap)
{
    const auto aTilesX = EpipTiles1D(aP0.x(),aP1.x(),aSz.x(),aOverlap.x());
    const auto aTilesY = EpipTiles1D(aP0.y(),aP1.y(),aSz.y(),aOverlap.y());
    std::vector<cRect2> aRes;
    for (const auto & aTY : aTilesY)
        for (const auto & aTX : aTilesX)
            aRes.push_back(cRect2(cPt2di(aTX.first,aTY.first),cPt2di(aTX.second,aTY.second)));
    return aRes;
}


// ============================================================
//  cEpipCropInfo : file describing a pair of crops
// ============================================================

void cEpipCropInfo::AddData(const cAuxAr2007 &anAux)
{
    MMVII::AddData(cAuxAr2007("ImageMaster", anAux), mNameImMaster);
    MMVII::AddData(cAuxAr2007("ImageSlave", anAux), mNameImSlave);
    MMVII::AddData(cAuxAr2007("CropMasterP0", anAux), mCropMaster0);
    MMVII::AddData(cAuxAr2007("CropMasterP1", anAux), mCropMaster1);
    MMVII::AddData(cAuxAr2007("CropSlaveP0", anAux), mCropSlave0);
    MMVII::AddData(cAuxAr2007("CropSlaveP1", anAux), mCropSlave1);
    MMVII::AddData(cAuxAr2007("SizeMaster", anAux), mSizeMaster);
    MMVII::AddData(cAuxAr2007("SizeSlave", anAux), mSizeSlave);
    MMVII::AddData(cAuxAr2007("FrameSizeMaster", anAux), mFrameSizeMaster);
    MMVII::AddData(cAuxAr2007("FrameSizeSlave", anAux), mFrameSizeSlave);
    MMVII::AddData(cAuxAr2007("Shift", anAux), mShift);
    MMVII::AddData(cAuxAr2007("ZInterval", anAux), mZInterval);
    MMVII::AddData(cAuxAr2007("DispRangeSlaveMinusMaster", anAux), mDispRange);
}

void AddData(const cAuxAr2007 &anAux, cEpipCropInfo &anInfo)
{
    anInfo.AddData(anAux);
}

void cEpipCropInfo::ToFile(const std::string &aNameFile) const
{
    // Same formatting precautions as the model : no truncation of small values
    PushPrecTxtSerial(std::numeric_limits<tREAL8>::max_digits10);
    SetFixedFloatTxtSerial(false);
    SaveInFile(*this, aNameFile);
    SetFixedFloatTxtSerial(true);
    PopPrecTxtSerial();
}

cEpipCropInfo cEpipCropInfo::FromFile(const std::string &aNameFile)
{
    cEpipCropInfo anInfo;
    ReadFromFile(anInfo, aNameFile);
    return anInfo;
}

cEpipCropInfo MakeEpipCropInfo(const cPt2di & aMasterP0, const cPt2di & aMasterP1,
                               const cEpipSlaveCrop & aSlave, const cPt2dr & aZIntv,
                               const cPt2di & aFrameSizeMaster, const cPt2di & aFrameSizeSlave,
                               const std::string & aNameImMaster, const std::string & aNameImSlave)
{
    cEpipCropInfo anInfo;
    anInfo.mNameImMaster = aNameImMaster;
    anInfo.mNameImSlave = aNameImSlave;
    anInfo.mCropMaster0 = aMasterP0;
    anInfo.mCropMaster1 = aMasterP1;
    anInfo.mCropSlave0 = aSlave.mP0;
    anInfo.mCropSlave1 = aSlave.mP1;
    anInfo.mSizeMaster = aMasterP1 - aMasterP0;
    anInfo.mSizeSlave = aSlave.mP1 - aSlave.mP0;
    anInfo.mFrameSizeMaster = aFrameSizeMaster;
    anInfo.mFrameSizeSlave = aFrameSizeSlave;
    anInfo.mShift = aSlave.mP0 - aMasterP0;
    anInfo.mZInterval = aZIntv;
    // x_slave - x_master (no crop) becomes the difference of the cropped coordinates
    const tREAL8 aOffset = aMasterP0.x() - aSlave.mP0.x();
    anInfo.mDispRange = cPt2dr(aSlave.mDispRange.x() + aOffset, aSlave.mDispRange.y() + aOffset);
    return anInfo;
}


// ============================================================
//  cEpipolarRectification
// ============================================================

cEpipolarRectification::cEpipolarRectification(const cSensorImage& aCam1,
                                               const cSensorImage& aCam2,
                                               const cParams&      aParams)
    : mCam1  (aCam1)
    , mCam2  (aCam2)
    , mParams(aParams)
{}

// ============================================================
//  Compute  (Algorithm 1 of the paper)
// ============================================================

cEpipPolyModel cEpipolarRectification::Compute()
{
    // ----------------------------------------------------------
    //  Step 1 – generate H-compatible pairs (both directions)
    //
    //  Forward  (master=1, slave=2) : gives center of I1 points
    //                                 and epipolar direction in I2
    //  Backward (master=2, slave=1) : gives center of I2 points
    //                                 and epipolar direction in I1
    // ----------------------------------------------------------

    std::vector<cEpiPair> aPairsATrain, aPairsATest, aPairsBTrain, aPairsBTest;
    cPt2dr aCenter1, aCenter2;
    cPt2dr aDir1,    aDir2;
    cPt2dr aZInterval1, aZInterval2;
    tREAL8 aGridStep1, aGridStep2;
    int aNbStepX1, aNbStepY1, aNbStepX2, aNbStepY2;

    GenerateData(mCam1, mCam2, aPairsATrain, aPairsATest, aCenter1, aDir2, aZInterval1, aGridStep1, aNbStepX1, aNbStepY1);
    GenerateData(mCam2, mCam1, aPairsBTrain, aPairsBTest, aCenter2, aDir1, aZInterval2, aGridStep2, aNbStepX2, aNbStepY2);

    // We must inverse aDir2 because it is computed in the direction from I2 to I1, but we want it in the direction from I1 to I2
    aDir2 = - aDir2;

    if ((aDir2.x() + aDir1.x()) <0)
    {
        aDir1 = -aDir1;
        aDir2 = -aDir2;
    }

    // Directions are already unit (GenerateData, which checks they are not null)
    aDir1 = VUnit(aDir1);
    aDir2 = VUnit(aDir2);
    // Training pairs only (excludes the held-out test pool).
    mNbPairs12 = aPairsATrain.size();
    mNbPairs21 = aPairsBTrain.size();

    // ----------------------------------------------------------
    //  Step 2 – apply rotation Rₖ to all points (eq. 25)
    //  q = (p - Ck) / Dk
    // ----------------------------------------------------------

    auto Rotate = [](const cPt2dr& p,
                     const cPt2dr& C,
                     const cPt2dr& D) -> cPt2dr
    {
        return (p - C) / D;
    };

    // aPairsA* already store (masterPt in I1, slavePt in I2)
    // aPairsB* store          (masterPt in I2, slavePt in I1)
    //   => swap to keep the convention (pt1, pt2)

    std::vector<cEpiPair> aRotPairsTrain, aRotPairsTest;
    aRotPairsTrain.reserve(aPairsATrain.size() + aPairsBTrain.size());
    aRotPairsTest.reserve(aPairsATest.size() + aPairsBTest.size());

    for (const auto& pr : aPairsATrain)
        aRotPairsTrain.push_back({ Rotate(pr.mP1, aCenter1, aDir1),
                                   Rotate(pr.mP2, aCenter2, aDir2) });
    for (const auto& pr : aPairsBTrain)
        aRotPairsTrain.push_back({ Rotate(pr.mP2, aCenter1, aDir1),   // I1 pt
                                   Rotate(pr.mP1, aCenter2, aDir2) }); // I2 pt

    for (const auto& pr : aPairsATest)
        aRotPairsTest.push_back({ Rotate(pr.mP1, aCenter1, aDir1),
                                  Rotate(pr.mP2, aCenter2, aDir2) });
    for (const auto& pr : aPairsBTest)
        aRotPairsTest.push_back({ Rotate(pr.mP2, aCenter1, aDir1),
                                  Rotate(pr.mP1, aCenter2, aDir2) });

    // ----------------------------------------------------------
    //  Step 3 – estimate V1 (with Y-axis identity) and V2 -- train pool only
    // ----------------------------------------------------------
    cPolyXY_N<tREAL8> aV1(mParams.mPolyDegree);
    cPolyXY_N<tREAL8> aV2(mParams.mPolyDegree);
    EstimateForwardPolynomials(aRotPairsTrain, aV1, aV2);

    // ----------------------------------------------------------
    //  Step 4 – estimate inverse polynomials W1, W2 -- train pool only
    // ----------------------------------------------------------

    cPolyXY_N<tREAL8> aW1(mParams.mPolyDegreeInv);
    cPolyXY_N<tREAL8> aW2(mParams.mPolyDegreeInv);
    EstimateInversePolynomial(aRotPairsTrain, aV1, aW1, UseFromPair::PT1);
    EstimateInversePolynomial(aRotPairsTrain, aV2, aW2, UseFromPair::PT2);

    // ----------------------------------------------------------
    //  Independent residuals on the held-out test pool
    // ----------------------------------------------------------
    EstimateIndepResiduals(aRotPairsTest, aV1, aV2, aW1, aW2);

    auto anEpipPolyModel = cEpipPolyModel {
        std::make_unique<cEpipPolyMapping>(aV1,aW1,aCenter1,aDir1,aZInterval1,aGridStep1,aNbStepX1,aNbStepY1),
        std::make_unique<cEpipPolyMapping>(aV2,aW2,aCenter2,aDir2,aZInterval2,aGridStep2,aNbStepX2,aNbStepY2),
    };
    anEpipPolyModel.ComputeCommonFraming(mCam1.PixelDomain().Box(),mCam2.PixelDomain().Box(),mParams.mEpipFrm, mParams.mMargin);

    if (! mParams.mNoWarnings)
    {
        // Diagnostic, warning only : mapping folding over the overlap
        const auto aFold1 = EpipFoldCount(anEpipPolyModel.EpipMap1(),mCam1,mCam2,aZInterval1);
        const auto aFold2 = EpipFoldCount(anEpipPolyModel.EpipMap2(),mCam2,mCam1,aZInterval2);
        for (const auto & [aNum,aFold] : {std::make_pair(1,aFold1),std::make_pair(2,aFold2)})
            if (aFold.first > 0)
                MMVII_USER_WARNING("The epipolar mapping of image " + ToStr(aNum) + " folds at " + ToStr(aFold.first) + " of "
                                   + ToStr(aFold.second) + " points tested in the overlap : the resampled image is unusable there, "
                                   + "try a lower Degree or check the Z interval");
    }

    return anEpipPolyModel;
}

std::pair<int,int> EpipFoldCount(const cEpipolarMapping & aMap, const cSensorImage & aSIM, const cSensorImage & aSIS, const cPt2dr & aZIntv)
{
    const tREAL8 aZMid = (aZIntv.x()+aZIntv.y()) / 2.0;
    const cPt2di aSz = aSIM.Sz();
    const int aNb = 20;
    int aNbPos = 0, aNbNeg = 0;
    for (int aKx=0 ; aKx<aNb ; aKx++)
        for (int aKy=0 ; aKy<aNb ; aKy++)
        {
            const cPt2dr aP((aKx+0.5)*aSz.x()/aNb,(aKy+0.5)*aSz.y()/aNb);
            if (! aSIS.IsVisibleOnImFrame(aSIS.Ground2Image(aSIM.ImageAndZ2Ground(TP3z(aP,aZMid)))))
                continue;
            const cPt2dr aQ = aMap.Value(aP);
            const cPt2dr aDx = aMap.Value(aP+cPt2dr(1,0)) - aQ;
            const cPt2dr aDy = aMap.Value(aP+cPt2dr(0,1)) - aQ;
            ((aDx.x()*aDy.y() - aDx.y()*aDy.x()) > 0 ? aNbPos : aNbNeg)++;
        }
    return {std::min(aNbPos,aNbNeg),aNbPos+aNbNeg};
}

// ============================================================
//  GenerateData  (Algorithm 2 of the paper)
// ============================================================



// ============================================================
//  EstimateForwardPolynomials  (Section 3.2.4 of the paper)
//
//  Unknowns :
//    xFree1[0..nFree1-1]  : free coefficients of V1  (a >= 1)
//    x2    [0..n2-1]      : all coefficients of V2
//
//  Observation for each pair (q1, q2) :
//    FreeBasis_V1(q1) * xFree1  -  Basis_V2(q2) * x2
//       = -LockedContrib_V1(q1)
//
//  The locked contribution of V1 is simply  q1.y()
//  because C[0,1]=1 and all other C[0,b]=0.
// ============================================================

// Tikhonov damping strength, relative to each coefficient's own scale (shared by V1,V2,W1,W2).
static constexpr tREAL8 aTikhoRelWeight = 0.1;

void cEpipolarRectification::EstimateForwardPolynomials(
        const std::vector<cEpiPair>& aPairs,
        cPolyXY_Nd&           aV1,
        cPolyXY_Nd&           aV2)
{
    // To calcul V1 params, use an auxiliary Polynom class
    //   which implements the property V(0,y) = y
    cPolyXY_N_IdentityOnYAxis<double>  aV1IdOnY(aV1.Degree());
    const int nFree1 = aV1IdOnY.NbFreeCoeffs();
    const int n2     = aV2.NbCoeffs();
    const int nTotal = nFree1 + n2;


    cLeasSqtAA<double> aSolver(nTotal);
    // Per-column sum of squared basis values : q1,q2 are raw pixel coordinates,
    // not normalized, so damping needs each column's own scale (below).
    std::vector<double> aSomSqBasis(nTotal, 0.0);

    for (const auto& pr : aPairs)
    {
        const cPt2dr& q1 = pr.mP1;
        const cPt2dr& q2 = pr.mP2;

        cDenseVect<double> aCoeff(nTotal);
        for (int k = 0; k < nTotal; ++k) aCoeff(k) = 0.0;

        // Free part of V1  (positive, indices 0..nFree1-1)
        {
            const cDenseVect<double> fb = aV1IdOnY.FreeBasisVector(q1);
            for (int k = 0; k < nFree1; ++k)
                aCoeff(k) = fb(k);
        }

        // V2 part (negative, indices nFree1..nTotal-1)
        {
            const cDenseVect<double> b2 = aV2.BasisVector(q2);
            for (int k = 0; k < n2; ++k)
                aCoeff(nFree1 + k) = -b2(k);
        }

        for (int k = 0; k < nTotal; ++k)
            aSomSqBasis[k] += aCoeff(k) * aCoeff(k);

        // RHS = -locked contribution of V1 at q1 = -q1.y()
        const double aRHS = -aV1IdOnY.LockedContribution(q1);

        aSolver.PublicAddObservation(1.0, aCoeff, aRHS);
    }

    // Tikhonov damping of the top-degree (a+b==Degree) coefficients of V1 (free
    // ones) and V2, weighted by each column's own basis-value scale.
    {
        const int d = aV1.Degree();
        int aFreeIdx = 0;
        for (int a=0; a<=d; a++)
        {
            for (int b=0; b<=d-a; b++)
            {
                if (cPolyXY_N_IdentityOnYAxis<double>::IsFreeCoeff(a,b))
                {
                    if (a+b == d)
                        aSolver.AddObsFixVar(aTikhoRelWeight * std::sqrt(aSomSqBasis[aFreeIdx]), aFreeIdx, 0.0);
                    ++aFreeIdx;
                }
            }
        }
        for (int a=0; a<=d; a++)
        {
            const int idx = nFree1 + aV2.Index(a, d-a);
            aSolver.AddObsFixVar(aTikhoRelWeight * std::sqrt(aSomSqBasis[idx]), idx, 0.0);
        }
    }

    const cDenseVect<double> aSol = aSolver.PublicSolve();
    mV1V2Var = aSolver.VarCurSol();

    // Restore V1 : locked coefficients are already set in the
    // constructor of cPolyXY_N_IdentityOnYAxis; just fill free ones.
    aV1IdOnY.SetFreeCoeffsFromSolution(aSol, 0);

    // Restore V2
    aV2.SetFromSolution(aSol, nFree1);

    // Set V1 from auxiliary class
    aV1=aV1IdOnY;
}


// ------------------------------------------------------------
//  EstimateInversePolynomial
//
//  Observation (eq. 34) :  Wk( qk.x ,  Vk(qk) ) = qk.y
// ------------------------------------------------------------

void cEpipolarRectification::EstimateInversePolynomial(
        const std::vector<cEpiPair>& aPairs,
        const cPolyXY_Nd&     aVk,
        cPolyXY_Nd&           aWk,
        UseFromPair                  aUsePt)
{
    const int nCoeff = aWk.NbCoeffs();
    cLeasSqtAA<double> aSolver(nCoeff);
    std::vector<double> aSomSqBasis(nCoeff, 0.0); // see EstimateForwardPolynomials

    for (const auto& pr : aPairs)
    {
        const cPt2dr& qk = aUsePt == UseFromPair::PT1 ? pr.mP1 : pr.mP2;

        // Epipolar coordinates of qk
        const double u = qk.x();
        const double v = aVk.Eval(qk);   // v = Vk(qk)

        // Observation : Wk(u, v) = qk.y
        const cDenseVect<double> aCoeff = aWk.BasisVector(u, v);
        const double             aRHS   = qk.y();

        for (int k = 0; k < nCoeff; ++k)
            aSomSqBasis[k] += aCoeff(k) * aCoeff(k);

        aSolver.PublicAddObservation(1.0, aCoeff, aRHS);
    }

    // Tikhonov damping of the top-degree coefficients (see EstimateForwardPolynomials).
    {
        const int d = aWk.Degree();
        for (int a=0; a<=d; a++)
        {
            const int idx = aWk.Index(a, d-a);
            aSolver.AddObsFixVar(aTikhoRelWeight * std::sqrt(aSomSqBasis[idx]), idx, 0.0);
        }
    }

    const cDenseVect<double> aSol = aSolver.PublicSolve();
    if (aUsePt == UseFromPair::PT1)
    {
        mW1Var = aSolver.VarCurSol();
    } else {
        mW2Var = aSolver.VarCurSol();
    }
    aWk.SetFromSolution(aSol);
}



// Sample points with the same pixel step in X and Y (square cells, not a fixed
// step count) ; edges always covered, last row/column may be a short remainder.
// At least aMinNbPts points are produced.
static std::vector<cPt2dr> EqualStepXYGrid(const cPt2di & aSz, int aMinNbPts,
                                            int & aOutNbStepX, int & aOutNbStepY, tREAL8 & aOutStep)
{
    const tREAL8 aW = aSz.x();
    const tREAL8 aH = aSz.y();

    aOutStep = std::sqrt((aW * aH) / std::max(1,aMinNbPts));
    for (;;)
    {
        aOutNbStepX = std::max(1,(int)std::ceil(aW / aOutStep - 1e-6));
        aOutNbStepY = std::max(1,(int)std::ceil(aH / aOutStep - 1e-6));
        if ((aOutNbStepX+1) * (aOutNbStepY+1) >= aMinNbPts)
            break;
        aOutStep *= 0.9;
    }

    std::vector<cPt2dr> aRes;
    aRes.reserve((aOutNbStepX+1) * (aOutNbStepY+1));
    for (int aKx=0; aKx<=aOutNbStepX; aKx++)
    {
        tREAL8 aX = std::min(aW, aKx * aOutStep);
        for (int aKy=0; aKy<=aOutNbStepY; aKy++)
            aRes.push_back(cPt2dr(aX, std::min(aH, aKy * aOutStep)));
    }
    return aRes;
}

void cEpipolarRectification::GenerateData(const cSensorImage &aCamM,
                                          const cSensorImage &aCamS,
                                          std::vector<cEpiPair> &aOutPairsTrain,
                                          std::vector<cEpiPair> &aOutPairsTest,
                                          cPt2dr &aOutCenterM,
                                          cPt2dr &aOutDirS,
                                          cPt2dr &aZInterval,
                                          tREAL8 &aOutGridStep, int &aOutNbStepX, int &aOutNbStepY
                                          ) const {
    aOutPairsTrain.clear();
    aOutPairsTest.clear();
    aOutCenterM = cPt2dr(0, 0);
    aOutDirS = cPt2dr(0, 0);

    // Altitude range for this master camera : mZIntv > tie-point-derived > native.
    aZInterval = EpipEffectiveZInterval(mParams,aCamM,mCam1,mCam2,mCachedHomolZIntv);
    const double Zmin = aZInterval.x();
    const double Zmax = aZInterval.y();

    // Altitude step : NbZLevels levels for train, NbZLevels-1 for test
    const int nZ = mParams.mNbZLevels * 2 - 1;

    // Z interval between 2 successive steps
    auto aStepZ = (Zmax - Zmin) / (nZ - 1);

    // Equal-size XY steps (pixels) : at least 30x the total unknowns of V1,V2,W1,W2
    // (measured minimum on real data ; 10x left outliers up to 20px).
    const int aNbCoeffsV2 = cPolyXY_Nd::NbCoeffsForDegree(mParams.mPolyDegree);
    const int aNbCoeffsV1 = aNbCoeffsV2 - (mParams.mPolyDegree + 1);
    const int aNbCoeffsW  = cPolyXY_Nd::NbCoeffsForDegree(mParams.mPolyDegreeInv);
    const int aMinNbPts = 30 * (aNbCoeffsV1 + aNbCoeffsV2 + 2 * aNbCoeffsW);
    std::vector<cPt2dr> aVPts = EqualStepXYGrid(aCamM.Sz(), aMinNbPts, aOutNbStepX, aOutNbStepY, aOutGridStep);

    int nCentroid = 0;
    size_t aNbSampled = 0;   // (point, Z) samples of the master image tried
    for (const auto &pM : aVPts) {
        // ---- Regular Z sweep : nZ levels. Alternates train/test by index ; Extra levels are random
        for (int aKZ = 0; aKZ < nZ + 2; aKZ++) {
            double Z0 = Zmin + aKZ * aStepZ;
            if (aKZ >= nZ) {
                Z0 = RandInInterval(Zmin, Zmax);
            }
            const cPt3dr aGround0 = aCamM.ImageAndZ2Ground(TP3z(pM, Z0));
            const cPt2dr pS0 = aCamS.Ground2Image(aGround0);
            ++aNbSampled;
            if (!aCamS.IsVisibleOnImFrame(pS0))
                continue;
            aOutCenterM = aOutCenterM + pM;
            ++nCentroid;
            if (aKZ % 2 == 0)
                aOutPairsTrain.push_back({pM, pS0});
            else
                aOutPairsTest.push_back({pM, pS0});

            const double Z1 = Z0 + aStepZ;
            if (Z1 > Zmax)
                continue;
            const cPt3dr aGround1 = aCamM.ImageAndZ2Ground(TP3z(pM, Z1));
            const cPt2dr pS1 = aCamS.Ground2Image(aGround1);
            if (!aCamS.IsVisibleOnImFrame(pS1))
                continue;
            cPt2dr aDelta = pS1 - pS0;
            if (SqN2(aDelta) > 1e-16) {
                aDelta = VUnit(aDelta);
                aOutDirS = aOutDirS + aDelta;
            }
        }
    }

    // Overlap : share of the sampled (point, Z) of the master image visible in the slave image
    const size_t aNbSeen = aOutPairsTrain.size() + aOutPairsTest.size();
    const std::string aOverlapMes = "Images " + aCamM.NameImage() + " and " + aCamS.NameImage() + " : " + ToStr((int)aNbSeen) + " of "
                                    + ToStr((int)aNbSampled) + " sampled points of the first (over the Z interval " + ToStr(aZInterval)
                                    + ") are visible in the second";
    MMVII_INTERNAL_ASSERT_User(
        (aOutPairsTrain.size() > mParams.mMinNbPairs) && (aOutPairsTest.size() > mParams.mMinNbPairs), eTyUEr::eUnClassedError,
        "The images do not overlap enough (" + ToStr((int)mParams.mMinNbPairs) + " points needed in each pool) : " + aOverlapMes
        + ". Check the images, the orientations and the Z interval (ZIntv)");
    if ((aNbSeen < 0.2 * aNbSampled) && !mParams.mNoWarnings)
        MMVII_USER_WARNING("Small overlap between the images : " + aOverlapMes);
    aOutCenterM = aOutCenterM * (1.0 / nCentroid);
    // Sum of unit displacements is null if none was kept or they cancel : cannot normalize
    MMVII_INTERNAL_ASSERT_User(SqN2(aOutDirS) > 1e-16, eTyUEr::eUnClassedError,
                               "Epipolar direction undefined (no usable Z displacement in the slave image)");
    aOutDirS = VUnit(aOutDirS);
}

// ============================================================
//  EstimateIndepResiduals : evaluate the fitted V1,V2,W1,W2 on the held-out
//  test pairs.
// ============================================================

void cEpipolarRectification::EstimateIndepResiduals(
        const std::vector<cEpiPair>& aPairsTest,
        const cPolyXY_Nd& aV1, const cPolyXY_Nd& aV2,
        const cPolyXY_Nd& aW1, const cPolyXY_Nd& aW2)
{
    double aSumV = 0.0;
    double aSumW1 = 0.0;
    double aSumW2 = 0.0;
    for (const auto& pr : aPairsTest)
    {
        const cPt2dr& q1 = pr.mP1;
        const cPt2dr& q2 = pr.mP2;

        const double v1 = aV1.Eval(q1);
        const double v2 = aV2.Eval(q2);
        const double eV = v1 - v2;
        aSumV += eV * eV;

        const double eW1 = aW1.Eval(cPt2dr(q1.x(), v1)) - q1.y();
        aSumW1 += eW1 * eW1;

        const double eW2 = aW2.Eval(cPt2dr(q2.x(), v2)) - q2.y();
        aSumW2 += eW2 * eW2;
    }

    const size_t aN = aPairsTest.size();
    mV1V2VarIndep = aN ? (aSumV  / aN) : 0.0;
    mW1VarIndep   = aN ? (aSumW1 / aN) : 0.0;
    mW2VarIndep   = aN ? (aSumW2 / aN) : 0.0;
}

// ============================================================
//  Z interval of a master camera : mZIntv > tie-point-derived > aCamM's own native
//  (EpipEffectiveZInterval below). Overriding a lower-priority source only warns.
// ============================================================

static cPt2dr ZIntervalFromHomolPts(const cEpipolarRectification::cParams & aParams,const cSensorImage & aCam1,const cSensorImage & aCam2,std::optional<cPt2dr> & aCache)
{
    if (aCache)
        return *aCache;

    int aNbKept = 0;
    tREAL8 aZmin = 0.0;
    tREAL8 aZmax = 0.0;
    for (const auto & aCple : aParams.mHomolPts->SetH())
    {
        const tREAL8 aRes = aCam1.PixResInterBundle(aCple, aCam2);
        if (aRes > aParams.mTiePMaxRes)
            continue;

        const tREAL8 aZ = aCam1.PInterBundle(aCple, aCam2).z();
        if (aNbKept == 0)
        {
            aZmin = aZmax = aZ;
        }
        else
        {
            aZmin = std::min(aZmin, aZ);
            aZmax = std::max(aZmax, aZ);
        }
        ++aNbKept;
    }

    const cPt2dr aSz = aCam1.PixelDomain().Box().Sz();
    const int aMinNb = std::max(aParams.mTiePMinNbFloor,
                                 (int)std::ceil(aParams.mTiePMinNbRatio * std::sqrt(aSz.x() * aSz.y())));
    MMVII_INTERNAL_ASSERT_User(aNbKept >= aMinNb, eTyUEr::eUnClassedError,
        "Not enough tie points after residual filtering to infer Z interval (" + ToStr(aNbKept)
        + " < " + ToStr(aMinNb) + "); provide ZIntv=[Zmin,Zmax], relax TiePMaxRes/TiePMinNb*, or add tie points");

    MMVII_INTERNAL_ASSERT_User((aZmax - aZmin) > 1e-6, eTyUEr::eUnClassedError,
        "Tie points give a degenerate (near-flat) Z interval [" + ToStr(aZmin) + "," + ToStr(aZmax)
        + "]; provide ZIntv=[Zmin,Zmax] explicitly for this scene");

    const tREAL8 aMargin = aParams.mZMargin * (aZmax - aZmin);
    aCache = cPt2dr(aZmin - aMargin, aZmax + aMargin);
    return *aCache;
}

cPt2dr EpipEffectiveZInterval(const cEpipolarRectification::cParams & aParams,const cSensorImage & aCamM,
                              const cSensorImage & aCam1,const cSensorImage & aCam2,std::optional<cPt2dr> & aCache)
{
    cPt2dr aResult(0,0);

    if (aParams.mZIntv)
    {
        if (aCamM.HasIntervalZ() && !aParams.mNoWarnings)
        {
            MMVII_USER_WARNING("Provided ZIntv overrides sensor's own Z validity interval");
        }
        if (aParams.mHomolPts && !aParams.mNoWarnings)
        {
            MMVII_USER_WARNING("Provided ZIntv overrides tie-point-derived Z validity interval");
        }
        aResult = *aParams.mZIntv;
    }
    else if (aParams.mHomolPts)
    {
        aResult = ZIntervalFromHomolPts(aParams,aCam1,aCam2,aCache);
        if (aCamM.HasIntervalZ() && !aParams.mNoWarnings)
        {
            MMVII_USER_WARNING("Tie-point-derived Z validity interval overrides sensor's own Z validity interval");
        }
    }
    else
    {
        MMVII_INTERNAL_ASSERT_User(aCamM.HasIntervalZ(), eTyUEr::eUnClassedError,
            "Sensor has no Z validity interval (no RPC); provide ZIntv=[Zmin,Zmax] or TieP=<dir>");
        aResult = aCamM.GetIntervalZ();
    }

    return aResult;
}

void cEpipolarModel::ComputeCommonFraming(
    const cTplBox<tREAL8,2> aBox1,
    const cTplBox<tREAL8,2> aBox2,
    eEpipFrm aFrmType,
    int aMargin
    )
{
    auto frame1 = EpipMap1().BoxOfFrontier(aBox1,1.0);
    auto frame2 = EpipMap2().BoxOfFrontier(aBox2,1.0);

    auto P1_0 = frame1.P0()-cPt2dr(aMargin,aMargin);
    auto P1_1 = frame1.P1()+cPt2dr(aMargin,aMargin);
    auto P2_0 = frame2.P0()-cPt2dr(aMargin,aMargin);
    auto P2_1 = frame2.P1()+cPt2dr(aMargin,aMargin);

    double yMin = std::max(P1_0.y(), P2_0.y());
    double yMax = std::min(P1_1.y(), P2_1.y());
    switch(aFrmType) {
    case eEpipFrm::eIntersect:
        yMin = std::max(P1_0.y(), P2_0.y());
        yMax = std::min(P1_1.y(), P2_1.y());
        break;
    case eEpipFrm::eUnion:
        yMin = std::min(P1_0.y(), P2_0.y());
        yMax = std::max(P1_1.y(), P2_1.y());
        break;
    case eEpipFrm::eImg_1:
        yMin = P1_0.y();
        yMax = P1_1.y();
        break;
    case eEpipFrm::eImg_2:
        yMin = P2_0.y();
        yMax = P2_1.y();
        break;
    case eEpipFrm::eNbVals:
    default:
        MMVII_INTERNAL_ERROR("Invalid value for FrameType : " + ToStr(aFrmType));
        break;
    }

    P1_0.y() = P2_0.y() = yMin;
    P1_1.y() = P2_1.y() = yMax;

    GetEpipMap1().SetEpipImFrame(cTplBox<double,2>(P1_0,P1_1).ToI());
    GetEpipMap2().SetEpipImFrame(cTplBox<double,2>(P2_0,P2_1).ToI());
}


void BenchEpipolar(cParamExeBench & aParam)
{
    if (! aParam.NewBench("Epipolar")) return;

    const std::string & aInDir = cMMVII_Appli::CurrentAppli().InputDirTestMMVII() + "/Epipolar/";
    const std::string & aTmpDir = cMMVII_Appli::CurrentAppli().TmpDirTestMMVII();

    const auto Name1 = std::string("Sensor1");
    const auto Name2 = std::string("Sensor2");

    auto aSensor1 = std::unique_ptr<cSensorImage>(ReadExternalSensor(aInDir + "RPC-" + Name1 + ".xml", Name1,false));
    auto aSensor2 = std::unique_ptr<cSensorImage>(ReadExternalSensor(aInDir + "RPC-" + Name2 + ".xml", Name2,false));

    // Epipolar geometry computing test
    auto aParams = cEpipolarRectification::cParams{5,9,3};
    auto aRectifier = cEpipolarRectification(*aSensor1, *aSensor2, aParams);
    auto aEpipModel = aRectifier.Compute();

    // The epipolar mappings do not fold over the overlap of the two images
    for (const auto& [aMap,aSM,aSS] : {std::make_tuple(&aEpipModel.EpipMap1(),aSensor1.get(),aSensor2.get()),
                                       std::make_tuple(&aEpipModel.EpipMap2(),aSensor2.get(),aSensor1.get())})
    {
        const auto aFold = EpipFoldCount(*aMap,*aSM,*aSS,aMap->ZInterval());
        MMVII_INTERNAL_ASSERT_bench((aFold.second>0) && (aFold.first==0), "Epipolar mapping folds (or no point tested)");
    }

    // Independent residuals : sanity only, finite and non-negative.
    MMVII_INTERNAL_ASSERT_bench(aRectifier.V1V2VarIndep() >= 0, "V1V2VarIndep is negative");
    MMVII_INTERNAL_ASSERT_bench(aRectifier.W1VarIndep() >= 0, "W1VarIndep is negative");
    MMVII_INTERNAL_ASSERT_bench(aRectifier.W2VarIndep() >= 0, "W2VarIndep is negative");
    MMVII_INTERNAL_ASSERT_bench(std::sqrt(aRectifier.V1V2VarIndep()) < 1.0, "V1V2VarIndep implausibly high");
    MMVII_INTERNAL_ASSERT_bench(std::sqrt(aRectifier.W1VarIndep()) < 1.0, "W1VarIndep implausibly high");
    MMVII_INTERNAL_ASSERT_bench(std::sqrt(aRectifier.W2VarIndep()) < 1.0, "W2VarIndep implausibly high");

    // Serialization round-trip of the pair model (as actually shipped).
    {
        const std::string aOriName = "OriTest";
        auto aModelFile = aTmpDir + "EpipModel." + GlobTaggedNameDefSerial();
        cEpipPairModel(aEpipModel.EpipMap1().Clone(), aEpipModel.EpipMap2().Clone(), aOriName, Name1, Name2).ToFile(aModelFile);
        auto aReloaded = cEpipPairModel::FromFile(aModelFile);

        MMVII_INTERNAL_ASSERT_bench(aReloaded.OriName() == aOriName, "Epip model round-trip : OriName mismatch");
        MMVII_INTERNAL_ASSERT_bench((aReloaded.ImName(1) == Name1) && (aReloaded.ImName(2) == Name2), "Epip model round-trip : image names mismatch");

        for (const auto& [aOrig,aReload] : {std::make_pair(static_cast<const cEpipPolyMapping*>(&aEpipModel.EpipMap1()),static_cast<const cEpipPolyMapping*>(&aReloaded.Map(1))),
                                             std::make_pair(static_cast<const cEpipPolyMapping*>(&aEpipModel.EpipMap2()),static_cast<const cEpipPolyMapping*>(&aReloaded.Map(2)))})
        {
            MMVII_INTERNAL_ASSERT_bench(aOrig->ZInterval() == aReload->ZInterval(), "Epip model round-trip : ZInterval mismatch");
            MMVII_INTERNAL_ASSERT_bench(aOrig->GridStep() == aReload->GridStep(), "Epip model round-trip : GridStep mismatch");
            MMVII_INTERNAL_ASSERT_bench(aOrig->NbStepX() == aReload->NbStepX(), "Epip model round-trip : NbStepX mismatch");
            MMVII_INTERNAL_ASSERT_bench(aOrig->NbStepY() == aReload->NbStepY(), "Epip model round-trip : NbStepY mismatch");
            MMVII_INTERNAL_ASSERT_bench(aOrig->EpipFrame().P0() == aReload->EpipFrame().P0(), "Epip model round-trip : EpipImFrame.P0 mismatch");
            MMVII_INTERNAL_ASSERT_bench(aOrig->EpipFrame().P1() == aReload->EpipFrame().P1(), "Epip model round-trip : EpipImFrame.P1 mismatch");
        }

        for (const auto& [aS,aOrig,aReload] : {std::make_tuple(&aSensor1,&aEpipModel.EpipMap1(),&aReloaded.Map(1)),
                                                std::make_tuple(&aSensor2,&aEpipModel.EpipMap2(),&aReloaded.Map(2))})
        {
            for (const auto& aPt : (*aS)->PtsSampledOnSensor(RandUnif_M_N(5,10),0))
            {
                auto aV0 = aOrig->Value(aPt);
                auto aV1 = aReload->Value(aPt);
                MMVII_INTERNAL_ASSERT_bench(Norm2(aV0-aV1) < 1e-8, "Epip model round-trip : Value mismatch after reload");
                auto aI0 = aOrig->Inverse(aV0);
                auto aI1 = aReload->Inverse(aV0);
                MMVII_INTERNAL_ASSERT_bench(Norm2(aI0-aI1) < 1e-8, "Epip model round-trip : Inverse mismatch after reload");
            }
        }
    }

    // RPCs in epipolar geometry computing test
    auto aEpipName1 = "Epip-" + Name1;
    auto aEpipName2 = "Epip-" + Name2;
    auto aEpipRPC1 = aTmpDir + aEpipName1 + ".xml";
    auto aEpipRPC2 = aTmpDir + aEpipName2 + ".xml";
    auto aResampSI1 = std::unique_ptr<cSensorImage>(aSensor1->GenerateSensorRPC(&aEpipModel.EpipMap1(), nullptr, false, aEpipName1));
    aResampSI1->ToFile(aEpipRPC1);
    auto aResampSI2 = std::unique_ptr<cSensorImage>(aSensor2->GenerateSensorRPC(&aEpipModel.EpipMap2(), nullptr, false, aEpipName2));
    aResampSI2->ToFile(aEpipRPC2);

    // Reread generated RPCs and check that they are consistent with the epipolar model
    auto aEpipSensor1 = std::unique_ptr<cSensorImage>(ReadExternalSensor(aEpipRPC1, aEpipName1,false));
    auto aEpipSensor2 = std::unique_ptr<cSensorImage>(ReadExternalSensor(aEpipRPC2, aEpipName2,false));
    aEpipSensor1->GetIntervalZ();

    // Test that epipolar image (mapping) of points is the same as the image of points by the epipolar RPC (direct then inverse)
    for (const auto& [aS,aES,aMap] : {std::make_tuple(&aSensor1, &aEpipSensor1, &aEpipModel.EpipMap1()),std::make_tuple(&aSensor2, &aEpipSensor2, &aEpipModel.EpipMap2())})
    {
        for (const auto& aPt : (*aS)->PtsSampledOnSensor(RandUnif_M_N(5,10),0))
        {
            auto aZ = RandInInterval((*aS)->GetIntervalZ());
            auto aPG = (*aS)->ImageAndZ2Ground(TP3z(ToR(aPt), aZ));
            auto aEpipPt = (*aES)->Ground2Image(aPG);
            auto aMapPt = aMap->Value(aPt);
            MMVII_INTERNAL_ASSERT_bench(Norm2(aMapPt - aEpipPt) < 0.1, "Epipolar geometry test failed (" + ToS(aMapPt) + " != " + ToS(aEpipPt));
        }
    }

    // Crops : same polynomials whatever the crop, so a ground point has one image position modulo the crop origin
    for (const auto& [aS,aMap] : {std::make_pair(&aSensor1, &aEpipModel.EpipMap1()),std::make_pair(&aSensor2, &aEpipModel.EpipMap2())})
    {
        const cPt2di aSzE = aMap->EpipImSz();
        const cPt2di aP0a = aMap->EpipFrame().P0() + cPt2di(aSzE.x()/4, aSzE.y()/5);
        const cPt2di aP0b = aMap->EpipFrame().P0() + cPt2di(aSzE.x()/2, aSzE.y()/3);
        const cPt2di aSzCrop = aSzE / 2;
        auto aFull = std::unique_ptr<cSensorImage>((*aS)->GenerateSensorRPC(aMap, nullptr, false, std::string("Crop-") + (*aS)->NameImage()));
        auto aCropA = std::unique_ptr<cSensorImage>(aFull->CropSensor(aP0a, aSzCrop));
        auto aCropB = std::unique_ptr<cSensorImage>(aFull->CropSensor(aP0b, aSzCrop));
        const std::string aNameA = aTmpDir + "CropA-" + (*aS)->NameImage() + ".xml";
        aCropA->ToFile(aNameA);
        auto aRereadA = std::unique_ptr<cSensorImage>(ReadExternalSensor(aNameA, "CropA-" + (*aS)->NameImage(), false));
        for (const auto& aPt : (*aS)->PtsSampledOnSensor(RandUnif_M_N(5,10),0))
        {
            auto aZ = RandInInterval(aMap->ZInterval());
            auto aPG = (*aS)->ImageAndZ2Ground(TP3z(ToR(aPt), aZ));
            // CropSensor is const : the source sensor keeps its own geometry
            MMVII_INTERNAL_ASSERT_bench(Norm2(aFull->Ground2Image(aPG) - aMap->Value(ToR(aPt))) < 0.1, "CropSensor modified its source sensor");
            auto aImA = aCropA->Ground2Image(aPG);
            auto aImB = aCropB->Ground2Image(aPG);
            // same ground point : positions differ exactly by the crop origins
            MMVII_INTERNAL_ASSERT_bench(Norm2((aImA + ToR(aP0a)) - (aImB + ToR(aP0b))) < 1e-6, "Crop RPC : offsets not the only difference (Ground2Image)");
            // and stay the epipolar position of the point, in the crop frame
            MMVII_INTERNAL_ASSERT_bench(Norm2((aImA + ToR(aP0a)) - aMap->Value(ToR(aPt))) < 0.1, "Crop RPC : Ground2Image mismatch");
            // analytic Jacobian path vs Ground2Image and vs finite differences : generated, cropped and reread RPC
            for (const cSensorImage * aSI : {aFull.get(), aCropA.get(), aRereadA.get()})
            {
                const auto aD = aSI->DiffGround2Im(aPG);
                const auto aF = aSI->DiffG2IByFiniteDiff(aPG);
                MMVII_INTERNAL_ASSERT_bench(Norm2(aD.mPIJ - aSI->Ground2Image(aPG)) < 1e-6, "RPC : DiffGround2Im value inconsistent with Ground2Image");
                const tREAL8 aNormGrad = Norm2(aD.mGradI) + Norm2(aD.mGradJ);
                MMVII_INTERNAL_ASSERT_bench(Norm2(aD.mGradI-aF.mGradI) + Norm2(aD.mGradJ-aF.mGradJ) < 1e-2 * aNormGrad, "RPC : DiffGround2Im gradient inconsistent with finite differences");
            }
            MMVII_INTERNAL_ASSERT_bench(Norm2(aRereadA->Ground2Image(aPG) - aImA) < 1e-3, "Crop RPC : file round trip changed Ground2Image");
            // pixel -> ground : same epipolar pixel in each crop's frame gives the same ground point
            auto aGA = aCropA->ImageAndZ2Ground(TP3z(aImA, aZ));
            auto aGB = aCropB->ImageAndZ2Ground(TP3z(aImB, aZ));
            MMVII_INTERNAL_ASSERT_bench(Norm2(aGA - aGB) < 1e-6 * (1.0 + Norm2(aPG)), "Crop RPC : ImageAndZ2Ground mismatch between crops");
            MMVII_INTERNAL_ASSERT_bench(Norm2(aGA - aPG) < 1e-3 * (1.0 + Norm2(aPG)), "Crop RPC : ImageAndZ2Ground does not give back the ground point");
        }
    }

    // Test that points with same y in master image have same y in slave image, for a random sampling of points on both sensors
    const std::pair<std::unique_ptr<cSensorImage>*, std::unique_ptr<cSensorImage>*> aPairs[] =
        {
            {&aEpipSensor1, &aEpipSensor2},
            {&aEpipSensor2, &aEpipSensor1}
        };
    for (const auto& [aES1,aES2] : aPairs)
    {
        for (const auto& aPt1 : (*aES1)->PtsSampledOnSensor(RandUnif_M_N(5,10),0))
        {
            auto aZ = RandInInterval((*aES1)->GetIntervalZ());
            auto aPG = (*aES1)->ImageAndZ2Ground(TP3z(ToR(aPt1), aZ));
            auto aPt2 = (*aES2)->Ground2Image(aPG);
            MMVII_INTERNAL_ASSERT_bench(std::abs(aPt2.y() - aPt1.y()) < 0.3, "Epipolar geometry test failed (y1:" + std::to_string(aPt1.y()) + " != y2:" + std::to_string(aPt2.y()));
        }
    }

    aParam.EndBench();
    return;
}


// ----------------------------------------------------------------------
//  BenchEpipolarResampling : crop + validity mask of ResampleEpipImage
//  (the code shared by EpipRectification and EpipResampling).
//  The source image is small (mask is only meaningful where the source ends), constant, sparse content.
// ----------------------------------------------------------------------

static void EpipolarResamplingBody(cParamExeBench & aParam)
{

    const std::string & aInDir = cMMVII_Appli::CurrentAppli().InputDirTestMMVII() + "/Epipolar/";
    const std::string & aTmpDir = cMMVII_Appli::CurrentAppli().TmpDirTestMMVII();

    const std::string aNameSens = "Sensor1";
    auto aSensor = std::unique_ptr<cSensorImage>(ReadExternalSensor(aInDir + "RPC-" + aNameSens + ".xml", aNameSens, false));
    auto aSensor2 = std::unique_ptr<cSensorImage>(ReadExternalSensor(aInDir + "RPC-Sensor2.xml", "Sensor2", false));
    auto aRectifier = cEpipolarRectification(*aSensor, *aSensor2, cEpipolarRectification::cParams{5,9,3});
    auto aEpipModel = aRectifier.Compute();
    const cEpipolarMapping & aMap = aEpipModel.EpipMap1();

    // Constant source image, much smaller than the sensor domain : most of the epipolar frame is outside it
    const cPt2di aSzSrc(1200,1000);
    const std::string aSrcName = aTmpDir + "EpipResSrc.tif";
    {
        cIm2D<tU_INT1> aSrc(aSzSrc);
        for (const auto & aPix : aSrc.DIm())
            aSrc.DIm().SetV(aPix,100);
        aSrc.DIm().ToFile(aSrcName);
    }

    // Crop = epipolar bounding box of the source frame, with a margin
    cPt2dr aMin(1e30,1e30), aMax(-1e30,-1e30);
    for (int aKx=0 ; aKx<=4 ; aKx++)
        for (int aKy=0 ; aKy<=4 ; aKy++)
        {
            cPt2dr aP = aMap.Value(cPt2dr(aSzSrc.x()*aKx/4.0, aSzSrc.y()*aKy/4.0));
            aMin = cPt2dr(std::min(aMin.x(),aP.x()),std::min(aMin.y(),aP.y()));
            aMax = cPt2dr(std::max(aMax.x(),aP.x()),std::max(aMax.y(),aP.y()));
        }
    const int aMarg = 20;
    cEpipCropMaskOpts aOpt;
    aOpt.mCropP0 = cPt2di(round_down(aMin.x())-aMarg, round_down(aMin.y())-aMarg);
    aOpt.mCropP1 = cPt2di(round_up(aMax.x())+aMarg, round_up(aMax.y())+aMarg);
    aOpt.mHasCrop = true;
    aOpt.mMaskOn = true;
    const cPt2di aSzCrop = aOpt.mCropP1 - aOpt.mCropP0;
    MMVII_INTERNAL_ASSERT_bench((aSzCrop.x()>0) && (aSzCrop.y()>0) && (aSzCrop.x()<6000) && (aSzCrop.y()<6000), "EpipolarResampling : unexpected crop size");

    std::unique_ptr<cInterpolator1D> anInterp(cDiffInterpolator1D::AllocFromNames({"Cubic","-0.5"}));
    std::vector<std::string> aToRemove {aSrcName};
    auto Path = [&](const std::string & aName) { aToRemove.push_back(aTmpDir+aName); return aTmpDir+aName; };

    // Run 1 : image + mask (default name pattern) + RPC
    const std::string aBase1 = "EpipRes1.tif";
    std::unique_ptr<cSensorImage> aEpipSI1(EpipSensor(*aSensor, aMap, aTmpDir + aBase1));
    ResampleEpipImage(aOpt, aSrcName, aMap, aEpipSI1.get(), *anInterp, false, aTmpDir, aBase1);
    const std::string aImName1 = Path(aBase1), aMaskName1 = Path("mask_EpipRes1.tif"), aRPCName1 = Path("RPC_EpipRes1.tif.xml");
    MMVII_INTERNAL_ASSERT_bench(ExistFile(aImName1) && ExistFile(aMaskName1) && ExistFile(aRPCName1), "EpipolarResampling : output files missing (default mask name)");

    auto aIm1 = cIm2D<tU_INT1>::FromFile(aImName1);
    auto aMask1 = cIm2D<tU_INT1>::FromFile(aMaskName1);
    MMVII_INTERNAL_ASSERT_bench((aIm1.DIm().Sz()==aSzCrop) && (aMask1.DIm().Sz()==aSzCrop), "EpipolarResampling : wrong output size");

    // Mask vs image : outside the mask the image is the default value, inside it is resampled data
    int aNbValid = 0;
    for (const auto & aPix : aMask1.DIm())
    {
        const int aM = aMask1.DIm().GetV(aPix);
        MMVII_INTERNAL_ASSERT_bench((aM==0) || (aM==1), "EpipolarResampling : mask is not binary");
        aNbValid += aM;
        MMVII_INTERNAL_ASSERT_bench((aIm1.DIm().GetV(aPix)==0) == (aM==0), "EpipolarResampling : image null <=> outside mask violated");
    }
    MMVII_INTERNAL_ASSERT_bench((aNbValid>0) && (aNbValid<aSzCrop.x()*aSzCrop.y()), "EpipolarResampling : mask should have valid and invalid pixels");

    // Mask vs geometry, through the sensors only : crop pixel -> ground (cropped RPC) -> source image
    auto aCropSI = std::unique_ptr<cSensorImage>(ReadExternalSensor(aRPCName1, "EpipRes1.tif", false));
    const tREAL8 aZ = (aMap.ZInterval().x() + aMap.ZInterval().y()) / 2.0;
    int aNbTested = 0;
    for (int aK=0 ; aK<300 ; aK++)
    {
        const cPt2di aPix(RandUnif_M_N(0,aSzCrop.x()-1), RandUnif_M_N(0,aSzCrop.y()-1));
        const cPt2dr aPSrc = aSensor->Ground2Image(aCropSI->ImageAndZ2Ground(TP3z(ToR(aPix),aZ)));
        const tREAL8 aMargin = std::min(std::min(aPSrc.x(), aSzSrc.x()-aPSrc.x()), std::min(aPSrc.y(), aSzSrc.y()-aPSrc.y()));
        if (std::abs(aMargin) < 3.0) continue;   // frontier : rounding may go either way
        aNbTested++;
        MMVII_INTERNAL_ASSERT_bench((aMargin>0) == (aMask1.DIm().GetV(aPix)==1), "EpipolarResampling : mask inconsistent with sensor geometry");
    }
    MMVII_INTERNAL_ASSERT_bench(aNbTested>100, "EpipolarResampling : too few pixels tested against geometry");

    // Run 2 : user mask pattern, no RPC : same mask, only the requested files
    aOpt.mMaskName = "m_$1_v.tif";
    ResampleEpipImage(aOpt, aSrcName, aMap, nullptr, *anInterp, false, aTmpDir, "EpipRes2.tif");
    const std::string aMaskName2 = Path("m_EpipRes2_v.tif");
    Path("EpipRes2.tif");
    MMVII_INTERNAL_ASSERT_bench(ExistFile(aMaskName2), "EpipolarResampling : mask with user pattern missing");
    MMVII_INTERNAL_ASSERT_bench(! ExistFile(aTmpDir + "RPC_EpipRes2.tif.xml"), "EpipolarResampling : RPC produced although not requested");
    auto aMask2 = cIm2D<tU_INT1>::FromFile(aMaskName2);
    for (const auto & aPix : aMask1.DIm())
        MMVII_INTERNAL_ASSERT_bench(aMask1.DIm().GetV(aPix)==aMask2.DIm().GetV(aPix), "EpipolarResampling : mask depends on the name pattern");

    // Run 3 : no mask requested : no mask file, same image
    aOpt.mMaskName = TheDefaultMaskNamePat;
    aOpt.mMaskOn = false;
    ResampleEpipImage(aOpt, aSrcName, aMap, nullptr, *anInterp, false, aTmpDir, "EpipRes3.tif");
    const std::string aImName3 = Path("EpipRes3.tif");
    MMVII_INTERNAL_ASSERT_bench(ExistFile(aImName3) && ! ExistFile(aTmpDir + "mask_EpipRes3.tif"), "EpipolarResampling : mask produced although not requested");
    auto aIm3 = cIm2D<tU_INT1>::FromFile(aImName3);
    for (const auto & aPix : aIm1.DIm())
        MMVII_INTERNAL_ASSERT_bench(aIm1.DIm().GetV(aPix)==aIm3.DIm().GetV(aPix), "EpipolarResampling : image changes with the mask option");

    // Run 4 : mask without image
    aOpt.mMaskOn = true;
    ResampleEpipImage(aOpt, aSrcName, aMap, nullptr, *anInterp, true, aTmpDir, "EpipRes4.tif");
    Path("mask_EpipRes4.tif");
    MMVII_INTERNAL_ASSERT_bench(ExistFile(aTmpDir + "mask_EpipRes4.tif") && ! ExistFile(aTmpDir + "EpipRes4.tif"), "EpipolarResampling : NoImage with mask");

    // Run 5 : the crop RPC made from the saved RPC of the full frame is the one fitted directly (same polynomials, bit for bit)
    cEpipCropMaskOpts aFullOpt;   // no crop, no mask
    std::unique_ptr<cSensorImage> aEpipSIFull(EpipSensor(*aSensor, aMap, aTmpDir + "EpipResFull.tif"));
    ResampleEpipImage(aFullOpt, aSrcName, aMap, aEpipSIFull.get(), *anInterp, true, aTmpDir, "EpipResFull.tif");
    const std::string aRPCFull = Path("RPC_EpipResFull.tif.xml");
    aOpt.mMaskOn = false;
    std::unique_ptr<cSensorImage> aEpipSI5(EpipSensor(*aSensor, aMap, aTmpDir + "EpipRes5.tif", aRPCFull));
    ResampleEpipImage(aOpt, aSrcName, aMap, aEpipSI5.get(), *anInterp, true, aTmpDir, "EpipRes5.tif");
    const std::string aRPCReused = Path("RPC_EpipRes5.tif.xml");
    auto aCropFit = std::unique_ptr<cSensorImage>(ReadExternalSensor(aRPCName1, "EpipRes1.tif", false));
    auto aCropReused = std::unique_ptr<cSensorImage>(ReadExternalSensor(aRPCReused, "EpipRes5.tif", false));
    for (int aK=0 ; aK<50 ; aK++)
    {
        const cPt3dr aG = aCropFit->ImageAndZ2Ground(TP3z(cPt2dr(RandInInterval(0.0,aSzCrop.x()-1.0),RandInInterval(0.0,aSzCrop.y()-1.0)),aZ));
        MMVII_INTERNAL_ASSERT_bench(Norm2(aCropFit->Ground2Image(aG) - aCropReused->Ground2Image(aG)) < 1e-9, "EpipolarResampling : RPC from the saved full-frame RPC differs from the fitted one");
    }

    // Generated RPC (full frame and crop) : integer image offsets and scales, as in NITF RPC00B
    for (const std::string & aNameRPC : {aRPCFull,aRPCReused})
    {
        std::ifstream aFile(aNameRPC);
        const std::string aText((std::istreambuf_iterator<char>(aFile)),std::istreambuf_iterator<char>());
        for (const std::string aTag : {"LINE_OFF","SAMP_OFF","LINE_SCALE","SAMP_SCALE"})
        {
            const size_t aPos = aText.find("<" + aTag + ">");
            MMVII_INTERNAL_ASSERT_bench(aPos!=std::string::npos, "EpipolarResampling : RPC tag " + aTag + " not found");
            const tREAL8 aVal = std::stod(aText.substr(aPos + aTag.size() + 2));
            MMVII_INTERNAL_ASSERT_bench(aVal==std::round(aVal), "EpipolarResampling : RPC " + aTag + " is not an integer");
        }
    }

    for (const auto & aName : aToRemove)
        RemoveFile(aName,SVP::Yes);

}

namespace {

// Catches MMVII_UserError/ASSERT_User via a throwing handler instead of abort();
// RAII-restored, same pattern as cProfileErrorCatcher.
struct cBenchNoZIntvError {};
std::string TheBenchLastErrorMes;   // message of the last error caught by cBenchErrorCatcher

void BenchNoZIntvErrorHandler(const std::string &, const std::string & aMes, const char *, int)
{
    TheBenchLastErrorMes = aMes;
    throw cBenchNoZIntvError{};
}

class cBenchErrorCatcher
{
public:
    cBenchErrorCatcher() : mPrev(MMVVI_Error) { MMVII_SetErrorHandler(BenchNoZIntvErrorHandler); }
    ~cBenchErrorCatcher() { MMVII_SetErrorHandler(mPrev); }
    cBenchErrorCatcher(const cBenchErrorCatcher &) = delete;
private:
    PtrMMVII_Error_Handler mPrev;
};

// Synthetic conic camera looking from aCenter to aTarget. Not cCamSimul: its
// terrestrial poses are too narrow a footprint for a well-conditioned RPC fit.
cSensorCamPC * BuildLookAtConicCam(const std::string & aName, const cPt3dr & aCenter,
                                    const cPt3dr & aTarget, cPerspCamIntrCalib * aCalib)
{
    const cPt3dr aK = VUnit(aTarget - aCenter);
    const cPt3dr aWorldUp(0,0,1);
    cPt3dr aI = VUnit(aK ^ ((std::abs(aK.z()) > 0.9) ? cPt3dr(1,0,0) : aWorldUp));
    cPt3dr aJ = VUnit(aK ^ aI);
    aI = aJ ^ aK;
    cRotation3D<tREAL8> aRot(M3x3FromCol(aI,aJ,aK),false);
    return new cSensorCamPC(aName,cIsometry3D<tREAL8>(aCenter,aRot),aCalib);
}

} // anonymous namespace


// ============================================================
//  BenchEpipolarSlaveCrop : the slave crop of a master crop must contain, over the whole Z
//  interval, every point matching a point of the master crop (same row, disparity in range).
// ============================================================

static void EpipolarSlaveCropBody(cParamExeBench & aParam)
{

    const std::string & aInDir = cMMVII_Appli::CurrentAppli().InputDirTestMMVII() + "/Epipolar/";
    auto aSensor1 = std::unique_ptr<cSensorImage>(ReadExternalSensor(aInDir + "RPC-Sensor1.xml", "Sensor1", false));
    auto aSensor2 = std::unique_ptr<cSensorImage>(ReadExternalSensor(aInDir + "RPC-Sensor2.xml", "Sensor2", false));
    auto aEpipModel = cEpipolarRectification(*aSensor1, *aSensor2, cEpipolarRectification::cParams{5,9,3}).Compute();
    const cEpipolarMapping * aMaps[2] = {&aEpipModel.EpipMap1(), &aEpipModel.EpipMap2()};
    const cSensorImage * aSIs[2] = {aSensor1.get(), aSensor2.get()};

    const cPt2dr aZIntv = EpipPairZInterval(*aMaps[0], *aMaps[1]);

    for (int aKM=0 ; aKM<2 ; aKM++)    // master = image 1, then image 2
    {
        const cEpipolarMapping & aMapM = *aMaps[aKM];
        const cEpipolarMapping & aMapS = *aMaps[1-aKM];
        const cSensorImage & aSIM = *aSIs[aKM];
        const cSensorImage & aSIS = *aSIs[1-aKM];

        // A crop in the middle of the master frame
        const cPt2di aSzM = aMapM.EpipImSz();
        const cPt2di aP0(aSzM.x()/2 - 150, aSzM.y()/2 - 100);
        const cPt2di aP1 = aP0 + cPt2di(300,200);
        const cEpipSlaveCrop aCrop = EpipSlaveCrop(aMapM,aMapS,aSIM,aSIS,aP0,aP1,aZIntv);

        MMVII_INTERNAL_ASSERT_bench((aCrop.mP0.y()==aP0.y()) && (aCrop.mP1.y()==aP1.y()), "EpipolarSlaveCrop : rows differ from the master crop");
        // Z really matters on this data : the disparity range is not a point
        MMVII_INTERNAL_ASSERT_bench(aCrop.mDispRange.y() - aCrop.mDispRange.x() > 1.0, "EpipolarSlaveCrop : disparity range too small");
        MMVII_INTERNAL_ASSERT_bench((aCrop.mMeanParallax > 1.0) && (aCrop.mMeanParallax <= aCrop.mDispRange.y() - aCrop.mDispRange.x() + 1e-9),
                                    "EpipolarSlaveCrop : mean parallax over Z must be significant and within the disparity range");
        MMVII_INTERNAL_ASSERT_bench(aCrop.mP1.x() - aCrop.mP0.x() > 300, "EpipolarSlaveCrop : slave crop narrower than the master one");

        // Any point of the master crop, at any Z of the interval, falls in the slave crop on the same row
        int aNbTested = 0;
        for (int aK=0 ; aK<2000 ; aK++)
        {
            const cPt2dr aQM(RandInInterval(aP0.x(),aP1.x()-1.0), RandInInterval(aP0.y(),aP1.y()-1.0));
            const cPt2dr aPM = aMapM.Inverse(aQM);
            if (aSIM.PixelDomain().Insideness(aPM) <= 0) continue;
            // Extremes of the interval are drawn on purpose
            const tREAL8 aZ = (aK%4==0) ? aZIntv.x() : ((aK%4==1) ? aZIntv.y() : RandInInterval(aZIntv.x(),aZIntv.y()));
            const cPt2dr aQS = aMapS.Value(aSIS.Ground2Image(aSIM.ImageAndZ2Ground(TP3z(aPM,aZ))));
            aNbTested++;
            MMVII_INTERNAL_ASSERT_bench(std::abs(aQS.y()-aQM.y()) < 0.01, "EpipolarSlaveCrop : rows of matching points differ");
            MMVII_INTERNAL_ASSERT_bench((aQS.x()>=aCrop.mP0.x()) && (aQS.x()<aCrop.mP1.x()), "EpipolarSlaveCrop : matching point outside the slave crop");
            const tREAL8 aDisp = aQS.x() - aQM.x();
            MMVII_INTERNAL_ASSERT_bench((aDisp>=aCrop.mDispRange.x()-0.5) && (aDisp<=aCrop.mDispRange.y()+0.5), "EpipolarSlaveCrop : disparity outside its range");
        }
        MMVII_INTERNAL_ASSERT_bench(aNbTested>1000, "EpipolarSlaveCrop : too few points tested");
    }

    // Crop info file : round trip, and the disparity convention in the cropped images' coordinates
    {
        const std::string aTmpDir = cMMVII_Appli::CurrentAppli().TmpDirTestMMVII();
        const cPt2di aSzM = aMaps[0]->EpipImSz();
        const cPt2di aP0(aSzM.x()/2 - 150, aSzM.y()/2 - 100);
        const cPt2di aP1 = aP0 + cPt2di(300,200);
        const cEpipSlaveCrop aCrop = EpipSlaveCrop(*aMaps[0],*aMaps[1],*aSIs[0],*aSIs[1],aP0,aP1,aZIntv);
        const cEpipCropInfo anInfo = MakeEpipCropInfo(aP0,aP1,aCrop,aZIntv,aMaps[0]->EpipImSz(),aMaps[1]->EpipImSz(),"CropM.tif","CropS.tif");
        const std::string aNameFile = aTmpDir + "EpipCropInfo." + GlobTaggedNameDefSerial();
        anInfo.ToFile(aNameFile);
        const cEpipCropInfo aReloaded = cEpipCropInfo::FromFile(aNameFile);
        RemoveFile(aNameFile,SVP::Yes);

        MMVII_INTERNAL_ASSERT_bench((aReloaded.mNameImMaster==anInfo.mNameImMaster) && (aReloaded.mNameImSlave==anInfo.mNameImSlave), "EpipolarSlaveCrop : info image names");
        MMVII_INTERNAL_ASSERT_bench((aReloaded.mCropMaster0==aP0) && (aReloaded.mCropMaster1==aP1), "EpipolarSlaveCrop : info master crop");
        MMVII_INTERNAL_ASSERT_bench((aReloaded.mCropSlave0==aCrop.mP0) && (aReloaded.mCropSlave1==aCrop.mP1), "EpipolarSlaveCrop : info slave crop");
        MMVII_INTERNAL_ASSERT_bench((aReloaded.mSizeMaster==aP1-aP0) && (aReloaded.mSizeSlave==aCrop.mP1-aCrop.mP0), "EpipolarSlaveCrop : info image sizes");
        MMVII_INTERNAL_ASSERT_bench((aReloaded.mFrameSizeMaster==aMaps[0]->EpipImSz()) && (aReloaded.mFrameSizeSlave==aMaps[1]->EpipImSz()), "EpipolarSlaveCrop : info frame sizes");
        MMVII_INTERNAL_ASSERT_bench(aReloaded.mShift==aCrop.mP0-aP0, "EpipolarSlaveCrop : info shift");
        MMVII_INTERNAL_ASSERT_bench(aReloaded.mZInterval==aZIntv, "EpipolarSlaveCrop : info Z interval");
        MMVII_INTERNAL_ASSERT_bench(aReloaded.mDispRange==anInfo.mDispRange, "EpipolarSlaveCrop : info disparity range");

        // Pixel of the master crop -> matching pixel of the slave crop : disparity (cropped coordinates) in the file's range
        int aNbTested = 0;
        for (int aK=0 ; aK<500 ; aK++)
        {
            const cPt2dr aQM(RandInInterval(aP0.x(),aP1.x()-1.0), RandInInterval(aP0.y(),aP1.y()-1.0));
            const cPt2dr aPM = aMaps[0]->Inverse(aQM);
            if (aSIs[0]->PixelDomain().Insideness(aPM) <= 0) continue;
            const tREAL8 aZ = RandInInterval(aZIntv.x(),aZIntv.y());
            const cPt2dr aQS = aMaps[1]->Value(aSIs[1]->Ground2Image(aSIs[0]->ImageAndZ2Ground(TP3z(aPM,aZ))));
            const cPt2dr aInCropM = aQM - ToR(aReloaded.mCropMaster0);
            const cPt2dr aInCropS = aQS - ToR(aReloaded.mCropSlave0);
            aNbTested++;
            const tREAL8 aDisp = aInCropS.x() - aInCropM.x();
            MMVII_INTERNAL_ASSERT_bench((aDisp>=aReloaded.mDispRange.x()-0.5) && (aDisp<=aReloaded.mDispRange.y()+0.5), "EpipolarSlaveCrop : disparity in cropped coordinates outside the file's range");
            MMVII_INTERNAL_ASSERT_bench((aInCropS.x()>=0) && (aInCropS.x()<aReloaded.mCropSlave1.x()-aReloaded.mCropSlave0.x()), "EpipolarSlaveCrop : match outside the slave crop image");
        }
        MMVII_INTERNAL_ASSERT_bench(aNbTested>200, "EpipolarSlaveCrop : too few points tested for the info file");
    }

    // Intersection of two Z intervals
    cEpipPolyMapping aMapA(cPolyXY_Nd(1),cPolyXY_Nd(1),cPt2dr(0,0),cPt2dr(1,0),cPt2dr(0,10),1.0,1,1);
    cEpipPolyMapping aMapB(cPolyXY_Nd(1),cPolyXY_Nd(1),cPt2dr(0,0),cPt2dr(1,0),cPt2dr(5,20),1.0,1,1);
    MMVII_INTERNAL_ASSERT_bench(EpipPairZInterval(aMapA,aMapB)==cPt2dr(5,10), "EpipolarSlaveCrop : Z intervals intersection");

}


// ============================================================
//  BenchEpipolarTiles : tiles of a pair == global resampling of the same region.
//  Resampling is pointwise, so every tile (images, masks) must be identical to the same window of the
//  global resampling, and the RPC of every tile must keep the pixel <-> ground relation.
// ============================================================

static void EpipolarTilesBody(cParamExeBench & aParam)
{

    // ---- 1. Tile grid (pure function)
    for (const int aLen : {1,99,100,101,250,1000})
        for (const int aSz : {10,100,300})
            for (const int aOv : {0,5,60})
            {
                if (aOv >= aSz) continue;
                const int aK0 = -7;
                const auto aTiles = EpipTiles1D(aK0,aK0+aLen,aSz,aOv);
                const std::string aMsg = "EpipolarTiles : grid len=" + ToStr(aLen) + " sz=" + ToStr(aSz) + " ov=" + ToStr(aOv);
                MMVII_INTERNAL_ASSERT_bench(! aTiles.empty(), aMsg);
                MMVII_INTERNAL_ASSERT_bench((aTiles.front().first==aK0) && (aTiles.back().second==aK0+aLen), aMsg + " : ends not covered");
                for (size_t aK=0 ; aK<aTiles.size() ; aK++)
                {
                    MMVII_INTERNAL_ASSERT_bench(aTiles[aK].second-aTiles[aK].first == std::min(aSz,aLen), aMsg + " : tile size");
                    if (aK>0)
                    {
                        MMVII_INTERNAL_ASSERT_bench(aTiles[aK].first > aTiles[aK-1].first, aMsg + " : tiles order");
                        MMVII_INTERNAL_ASSERT_bench(aTiles[aK-1].second - aTiles[aK].first >= aOv, aMsg + " : overlap");
                    }
                }
                if (aLen > aSz)   // minimal number of tiles
                {
                    const int aStep = aSz - aOv;
                    MMVII_INTERNAL_ASSERT_bench((int)aTiles.size() == (aLen-aOv+aStep-1)/aStep, aMsg + " : number of tiles");
                }
            }
    {
        const cPt2di aP0(3,-5), aP1(263,195);
        const auto aTiles = EpipTiles(aP0,aP1,cPt2di(100,80),cPt2di(20,10));
        cIm2D<tU_INT1> aCover(aP1-aP0,nullptr,eModeInitImage::eMIA_Null);
        for (const auto & aTile : aTiles)
        {
            MMVII_INTERNAL_ASSERT_bench((aTile.P0().x()>=aP0.x()) && (aTile.P0().y()>=aP0.y()) && (aTile.P1().x()<=aP1.x()) && (aTile.P1().y()<=aP1.y()),
                                        "EpipolarTiles : 2D tile outside the region");
            for (const auto & aPix : aTile)
                aCover.DIm().SetV(aPix-aP0,1);
        }
        for (const auto & aPix : aCover.DIm())
            MMVII_INTERNAL_ASSERT_bench(aCover.DIm().GetV(aPix)==1, "EpipolarTiles : 2D region not covered");
        MMVII_INTERNAL_ASSERT_bench(aTiles.size()==3*3, "EpipolarTiles : 2D number of tiles");
    }

    // Overlap per axis : the 2D grid is the product of the two 1D grids
    {
        const cPt2di aP0(0,0), aP1(500,300), aSz(100,100);
        const cPt2di aOv(30,0);
        MMVII_INTERNAL_ASSERT_bench(EpipTiles(aP0,aP1,aSz,aOv).size() == EpipTiles1D(0,500,100,30).size()*EpipTiles1D(0,300,100,0).size(),
                                    "EpipolarTiles : 2D grid with an overlap per axis");
    }

    // Tile names : fixed width indices, inserted before the extension
    MMVII_INTERNAL_ASSERT_bench(EpipTileName("Epip_%1_%2.tif",1,2,3,3)=="Epip_%1_%2_t01_02.tif", "EpipolarTiles : tile name");
    MMVII_INTERNAL_ASSERT_bench((EpipNameWithExtension("Epip_%1_%2.tif")=="Epip_%1_%2.tif") && (EpipNameWithExtension("a.b.tiff")=="a.b.tiff"), "EpipolarTiles : name with extension changed");
    MMVII_INTERNAL_ASSERT_bench((EpipNameWithExtension("Epip_%1_%2")=="Epip_%1_%2.tif") && (EpipNameWithExtension("out.v2_%2")=="out.v2_%2.tif")
                                && (EpipNameWithExtension("m_$1")=="m_$1.tif") && (EpipNameWithExtension("name.")=="name.tif"), "EpipolarTiles : .tif not added to a name without extension");
    MMVII_INTERNAL_ASSERT_bench(EpipTileName("Out",0,5,2,6)=="Out_t00_05", "EpipolarTiles : tile name without extension");
    MMVII_INTERNAL_ASSERT_bench(EpipTileName("A.B.tif",7,15,120,40)=="A.B_t007_015.tif", "EpipolarTiles : tile name, wide indices");

    // Tile index file : round trip
    {
        cEpipTilesInfo anInfo;
        anInfo.mNameModel = "Model.xml";
        anInfo.mMaster = 2;
        anInfo.mSzTiles = cPt2di(600,500);
        anInfo.mSzOverL = cPt2di(60,40);
        anInfo.mRegion0 = cPt2di(-3,4);
        anInfo.mRegion1 = cPt2di(1300,1400);
        anInfo.mDispRange = cPt2dr(-12.25,1.5e-17);
        for (int aK=0 ; aK<3 ; aK++)
        {
            cEpipTileEntry anEntry;
            anEntry.mRow = aK/2;
            anEntry.mCol = aK%2;
            anEntry.mInfoFile = "Tile_" + ToStr(aK) + ".Info.xml";
            anInfo.mTiles.push_back(anEntry);
        }
        const std::string aNameFile = cMMVII_Appli::CurrentAppli().TmpDirTestMMVII() + "EpipTilesInfo." + GlobTaggedNameDefSerial();
        anInfo.ToFile(aNameFile);
        const cEpipTilesInfo aReloaded = cEpipTilesInfo::FromFile(aNameFile);
        RemoveFile(aNameFile,SVP::Yes);
        MMVII_INTERNAL_ASSERT_bench((aReloaded.mNameModel==anInfo.mNameModel) && (aReloaded.mMaster==2), "EpipolarTiles : index model or master");
        MMVII_INTERNAL_ASSERT_bench((aReloaded.mSzTiles==anInfo.mSzTiles) && (aReloaded.mSzOverL==anInfo.mSzOverL), "EpipolarTiles : index tile sizes");
        MMVII_INTERNAL_ASSERT_bench((aReloaded.mRegion0==anInfo.mRegion0) && (aReloaded.mRegion1==anInfo.mRegion1), "EpipolarTiles : index region");
        MMVII_INTERNAL_ASSERT_bench(aReloaded.mDispRange==anInfo.mDispRange, "EpipolarTiles : index disparity range");
        MMVII_INTERNAL_ASSERT_bench(aReloaded.mTiles.size()==3, "EpipolarTiles : index number of tiles");
        for (size_t aK=0 ; aK<3 ; aK++)
            MMVII_INTERNAL_ASSERT_bench((aReloaded.mTiles[aK].mRow==anInfo.mTiles[aK].mRow) && (aReloaded.mTiles[aK].mCol==anInfo.mTiles[aK].mCol)
                                        && (aReloaded.mTiles[aK].mInfoFile==anInfo.mTiles[aK].mInfoFile), "EpipolarTiles : index entry");
    }

    // ---- 2. Tiles of a pair against the global resampling
    const std::string & aInDir = cMMVII_Appli::CurrentAppli().InputDirTestMMVII() + "/Epipolar/";
    const std::string & aTmpDir = cMMVII_Appli::CurrentAppli().TmpDirTestMMVII();
    auto aSensor1 = std::unique_ptr<cSensorImage>(ReadExternalSensor(aInDir + "RPC-Sensor1.xml", "Sensor1", false));
    auto aSensor2 = std::unique_ptr<cSensorImage>(ReadExternalSensor(aInDir + "RPC-Sensor2.xml", "Sensor2", false));
    auto aEpipModel = cEpipolarRectification(*aSensor1, *aSensor2, cEpipolarRectification::cParams{5,9,3}).Compute();
    const cEpipolarMapping & aMap1 = aEpipModel.EpipMap1();
    const cEpipolarMapping & aMap2 = aEpipModel.EpipMap2();
    const cPt2dr aZIntv = EpipPairZInterval(aMap1,aMap2);

    // Source images : white noise (shifts, rounding and cubic overshoot all show), small but at the origin of the sensors.
    // The second one covers where the first one is seen in image 2.
    const cPt2di aSzSrc1(1200,1000);
    cPt2dr aMaxS(0,0);
    for (int aKx=0 ; aKx<=4 ; aKx++)
        for (int aKy=0 ; aKy<=4 ; aKy++)
            for (const tREAL8 aZ : {aZIntv.x(),aZIntv.y()})
            {
                const cPt2dr aP2 = aSensor2->Ground2Image(aSensor1->ImageAndZ2Ground(TP3z(cPt2dr(aSzSrc1.x()*aKx/4.0,aSzSrc1.y()*aKy/4.0),aZ)));
                aMaxS = cPt2dr(std::max(aMaxS.x(),aP2.x()),std::max(aMaxS.y(),aP2.y()));
            }
    const cPt2di aSzSrc2(round_up(aMaxS.x())+20,round_up(aMaxS.y())+20);
    MMVII_INTERNAL_ASSERT_bench((aSzSrc2.x()>100) && (aSzSrc2.y()>100) && (aSzSrc2.x()<5000) && (aSzSrc2.y()<5000), "EpipolarTiles : unexpected second source size");

    std::vector<std::string> aToRemove;
    auto MakeSrc = [&](const cPt2di & aSz,const std::string & aName,tU_INT4 aSeed)
    {
        cIm2D<tU_INT1> anIm(aSz);
        for (const auto & aPix : anIm.DIm())
        {
            tU_INT4 aH = tU_INT4(aPix.x())*2654435761u ^ tU_INT4(aPix.y())*40503u ^ aSeed;
            aH ^= aH>>13;
            aH *= 0x5bd1e995u;
            aH ^= aH>>15;
            anIm.DIm().SetV(aPix,int(aH&255));
        }
        anIm.DIm().ToFile(aName);
        aToRemove.push_back(aName);
    };
    const std::string aSrcName1 = aTmpDir + "EpipTileSrc1.tif", aSrcName2 = aTmpDir + "EpipTileSrc2.tif";
    MakeSrc(aSzSrc1,aSrcName1,12345u);
    MakeSrc(aSzSrc2,aSrcName2,67890u);

    // Region of the master epipolar frame : bounding box of the master source frame, with a margin
    cPt2dr aMin(1e30,1e30), aMax(-1e30,-1e30);
    for (int aKx=0 ; aKx<=4 ; aKx++)
        for (int aKy=0 ; aKy<=4 ; aKy++)
        {
            const cPt2dr aP = aMap1.Value(cPt2dr(aSzSrc1.x()*aKx/4.0, aSzSrc1.y()*aKy/4.0));
            aMin = cPt2dr(std::min(aMin.x(),aP.x()),std::min(aMin.y(),aP.y()));
            aMax = cPt2dr(std::max(aMax.x(),aP.x()),std::max(aMax.y(),aP.y()));
        }
    const int aMarg = 20;
    const cPt2di aR0(round_down(aMin.x())-aMarg, round_down(aMin.y())-aMarg);
    const cPt2di aR1(round_up(aMax.x())+aMarg, round_up(aMax.y())+aMarg);
    const auto aTiles = EpipTiles(aR0,aR1,cPt2di(600,600),cPt2di(60,60));
    MMVII_INTERNAL_ASSERT_bench((aTiles.size()>=4) && (aTiles.size()<=16), "EpipolarTiles : unexpected number of tiles " + ToStr((int)aTiles.size()) + ", region " + ToStr(aR0) + " " + ToStr(aR1));

    std::unique_ptr<cInterpolator1D> anInterp(cDiffInterpolator1D::AllocFromNames({"Cubic","-0.5"}));
    auto Resample = [&](const std::string & aSrc,const cEpipolarMapping & aMap,const cSensorImage * aSI,const cPt2di & aP0,const cPt2di & aP1,const std::string & aBase)
    {
        cEpipCropMaskOpts anOpt;
        anOpt.mHasCrop = true;
        anOpt.mMaskOn = true;
        anOpt.mCropP0 = aP0;
        anOpt.mCropP1 = aP1;
        std::unique_ptr<cSensorImage> aEpipSI(aSI ? EpipSensor(*aSI,aMap,aTmpDir+aBase) : nullptr);
        ResampleEpipImage(anOpt,aSrc,aMap,aEpipSI.get(),*anInterp,false,aTmpDir,aBase);
        aToRemove.push_back(aTmpDir+aBase);
        aToRemove.push_back(aTmpDir+"mask_"+LastPrefix(aBase)+".tif");
        aToRemove.push_back(aTmpDir+"RPC_"+aBase+".xml");
    };
    auto ReadIm = [&](const std::string & aName) { return cIm2D<tU_INT1>::FromFile(aTmpDir+aName); };
    auto ReadRPC = [&](const std::string & aBase) { return std::unique_ptr<cSensorImage>(ReadExternalSensor(aTmpDir+"RPC_"+aBase+".xml",aBase,false)); };
    // Tile image == window of the global one starting at aOff
    auto CmpWindow = [&](const cIm2D<tU_INT1> & aTile,const cIm2D<tU_INT1> & aGlob,const cPt2di & aOff,const std::string & aMsg)
    {
        const cPt2di aSzT = aTile.DIm().Sz(), aSzG = aGlob.DIm().Sz();
        MMVII_INTERNAL_ASSERT_bench((aOff.x()>=0) && (aOff.y()>=0) && (aOff.x()+aSzT.x()<=aSzG.x()) && (aOff.y()+aSzT.y()<=aSzG.y()), aMsg + " : window outside the global image");
        int aNbDif = 0;
        for (const auto & aPix : aTile.DIm())
            if (aTile.DIm().GetV(aPix) != aGlob.DIm().GetV(aPix+aOff))
                aNbDif++;
        MMVII_INTERNAL_ASSERT_bench(aNbDif==0, aMsg + " : " + ToStr(aNbDif) + " pixels differ from the global resampling");
    };

    // Reference reading the whole source (the resampling reads only the window it needs)
    auto CmpFullRead = [&](const std::string & aSrc,const cEpipolarMapping & aMap,const cPt2di & aP0,const cPt2di & aP1,const cIm2D<tU_INT1> & aTileIm,const std::string & aMsg)
    {
        MMVII_INTERNAL_ASSERT_bench(aTileIm.DIm().Sz()==aP1-aP0, aMsg + " : wrong size");
        const cEpipCropMapping aCropMap(aMap,ToR(aP0));
        std::unique_ptr<cDataGenUnTypedIm<2>> aSrcIm(ReadIm2DGen(aSrc));
        std::unique_ptr<cDataGenUnTypedIm<2>> aRef(aSrcIm->AllocReSampleGen(*anInterp,aCropMap,cRect2(cPt2di(0,0),aP1-aP0)));
        int aNbDif = 0;
        for (const auto & aPix : aTileIm.DIm())
            if (int(aRef->VD_GetV(aPix)) != aTileIm.DIm().GetV(aPix))
                aNbDif++;
        MMVII_INTERNAL_ASSERT_bench(aNbDif==0, aMsg + " : " + ToStr(aNbDif) + " pixels differ from the resampling of the whole source");
    };

    // Crops of the second image : derived from each master tile
    std::vector<cEpipSlaveCrop> aSlaveCrops;
    cPt2di aS0(1000000,1000000), aS1(-1000000,-1000000);
    for (const auto & aTile : aTiles)
    {
        aSlaveCrops.push_back(EpipSlaveCrop(aMap1,aMap2,*aSensor1,*aSensor2,aTile.P0(),aTile.P1(),aZIntv));
        const auto & aSC = aSlaveCrops.back();
        aS0 = cPt2di(std::min(aS0.x(),aSC.mP0.x()),std::min(aS0.y(),aSC.mP0.y()));
        aS1 = cPt2di(std::max(aS1.x(),aSC.mP1.x()),std::max(aS1.y(),aSC.mP1.y()));
    }

    // Global references : the whole region (master) and the union of the derived crops (slave)
    Resample(aSrcName1,aMap1,aSensor1.get(),aR0,aR1,"TileGlob1.tif");
    Resample(aSrcName2,aMap2,aSensor2.get(),aS0,aS1,"TileGlob2.tif");
    const auto aGlobIm1 = ReadIm("TileGlob1.tif"), aGlobMask1 = ReadIm("mask_TileGlob1.tif");
    const auto aGlobIm2 = ReadIm("TileGlob2.tif"), aGlobMask2 = ReadIm("mask_TileGlob2.tif");
    const auto aGlobRPC1 = ReadRPC("TileGlob1.tif");
    const auto aGlobRPC2 = ReadRPC("TileGlob2.tif");

    cIm2D<tU_INT1> aCanvas(aR1-aR0,nullptr,eModeInitImage::eMIA_Null);
    cIm2D<tU_INT1> aPainted(aR1-aR0,nullptr,eModeInitImage::eMIA_Null);
    for (size_t aKT=0 ; aKT<aTiles.size() ; aKT++)
    {
        const cRect2 & aTile = aTiles[aKT];
        const cEpipSlaveCrop & aSC = aSlaveCrops[aKT];
        const std::string aK = ToStr((int)aKT);
        const std::string aBaseM = "TileM_" + aK + ".tif", aBaseS = "TileS_" + aK + ".tif";
        Resample(aSrcName1,aMap1,aSensor1.get(),aTile.P0(),aTile.P1(),aBaseM);
        Resample(aSrcName2,aMap2,aSensor2.get(),aSC.mP0,aSC.mP1,aBaseS);

        // Images and masks : identical to the windows of the global resampling
        const auto aImM = ReadIm(aBaseM);
        const auto aImS = ReadIm(aBaseS);
        CmpWindow(aImM,aGlobIm1,aTile.P0()-aR0,"EpipolarTiles : master image of tile " + aK);
        CmpWindow(ReadIm("mask_TileM_"+aK+".tif"),aGlobMask1,aTile.P0()-aR0,"EpipolarTiles : master mask of tile " + aK);
        CmpWindow(aImS,aGlobIm2,aSC.mP0-aS0,"EpipolarTiles : slave image of tile " + aK);
        CmpWindow(ReadIm("mask_TileS_"+aK+".tif"),aGlobMask2,aSC.mP0-aS0,"EpipolarTiles : slave mask of tile " + aK);

        // First and last tiles : same as reading the whole source
        if ((aKT==0) || (aKT+1==aTiles.size()))
        {
            CmpFullRead(aSrcName1,aMap1,aTile.P0(),aTile.P1(),aImM,"EpipolarTiles : master tile " + aK);
            CmpFullRead(aSrcName2,aMap2,aSC.mP0,aSC.mP1,aImS,"EpipolarTiles : slave tile " + aK);
        }

        // Assembly of the master tiles : overlaps must agree exactly
        for (const auto & aPix : aImM.DIm())
        {
            const cPt2di aPC = aPix + aTile.P0() - aR0;
            if (aPainted.DIm().GetV(aPC))
                MMVII_INTERNAL_ASSERT_bench(aCanvas.DIm().GetV(aPC)==aImM.DIm().GetV(aPix), "EpipolarTiles : tiles disagree in their overlap");
            aCanvas.DIm().SetV(aPC,aImM.DIm().GetV(aPix));
            aPainted.DIm().SetV(aPC,1);
        }

        // RPC of the tiles : same ground relation as the global RPC and as the original sensor
        const auto aRPCM = ReadRPC(aBaseM);
        const auto aRPCS = ReadRPC(aBaseS);
        const cPt2di aOffM = aTile.P0() - aR0, aOffS = aSC.mP0 - aS0;
        const tREAL8 aTol = 0.01;   // px, order of the RPC fit accuracy (tile vs global RPC : same fit, 1e-6)
        for (int aKP=0 ; aKP<30 ; aKP++)
        {
            const cPt2dr aPCrop(RandInInterval(0.0,aTile.Sz().x()-1.0),RandInInterval(0.0,aTile.Sz().y()-1.0));
            const tREAL8 aZ = RandInInterval(aZIntv.x(),aZIntv.y());
            // inverse : tile RPC and global RPC give the same ground line
            const cPt3dr aGTile = aRPCM->ImageAndZ2Ground(TP3z(aPCrop,aZ));
            const cPt3dr aGGlob = aGlobRPC1->ImageAndZ2Ground(TP3z(aPCrop+ToR(aOffM),aZ));
            MMVII_INTERNAL_ASSERT_bench(Norm2(aSensor1->Ground2Image(aGTile)-aSensor1->Ground2Image(aGGlob)) < aTol, "EpipolarTiles : tile and global RPC inverse differ (master)");
            // direct : same image position up to the crop offset
            MMVII_INTERNAL_ASSERT_bench(Norm2(aRPCM->Ground2Image(aGTile) - (aGlobRPC1->Ground2Image(aGTile)-ToR(aOffM))) < 1e-6, "EpipolarTiles : tile and global RPC direct differ (master)");
            MMVII_INTERNAL_ASSERT_bench(Norm2(aRPCM->Ground2Image(aGTile) - aPCrop) < aTol, "EpipolarTiles : RPC direct o inverse is not the identity (master)");
            // relation with the original sensor : the epipolar position of the source pixel
            const cPt2dr aPSrc = aSensor1->Ground2Image(aGTile);
            MMVII_INTERNAL_ASSERT_bench(Norm2(aMap1.Value(aPSrc) - ToR(aTile.P0()) - aPCrop) < aTol, "EpipolarTiles : RPC not consistent with the epipolar mapping (master)");

            // slave side : the same ground point, seen by the slave tile RPC and by the global slave RPC
            const cPt2dr aPS = aRPCS->Ground2Image(aGTile);
            MMVII_INTERNAL_ASSERT_bench(Norm2(aPS - (aGlobRPC2->Ground2Image(aGTile)-ToR(aOffS))) < 1e-6, "EpipolarTiles : tile and global RPC direct differ (slave)");
            MMVII_INTERNAL_ASSERT_bench(Norm2(aMap2.Value(aSensor2->Ground2Image(aGTile)) - ToR(aSC.mP0) - aPS) < aTol, "EpipolarTiles : RPC not consistent with the epipolar mapping (slave)");

            // pair : same row in the two crops, disparity (cropped coordinates) inside the range of the info file
            const cEpipCropInfo anInfo = MakeEpipCropInfo(aTile.P0(),aTile.P1(),aSC,aZIntv,aMap1.EpipImSz(),aMap2.EpipImSz(),aBaseM,aBaseS);
            MMVII_INTERNAL_ASSERT_bench(std::abs(aPS.y()-aPCrop.y()) < aTol, "EpipolarTiles : matching points on different rows");
            const tREAL8 aDisp = aPS.x() - aPCrop.x();
            MMVII_INTERNAL_ASSERT_bench((aDisp>=anInfo.mDispRange.x()-0.5) && (aDisp<=anInfo.mDispRange.y()+0.5), "EpipolarTiles : disparity outside its range");
        }
    }
    // A crop entirely outside the source : null image, as when reading the whole source
    {
        const cPt2di aO0 = aR1 + cPt2di(200,200), aO1 = aO0 + cPt2di(150,100);
        Resample(aSrcName1,aMap1,nullptr,aO0,aO1,"TileOut.tif");
        const auto aImOut = ReadIm("TileOut.tif");
        for (const auto & aPix : aImOut.DIm())
            MMVII_INTERNAL_ASSERT_bench(aImOut.DIm().GetV(aPix)==0, "EpipolarTiles : crop outside the source is not null");
        CmpFullRead(aSrcName1,aMap1,aO0,aO1,aImOut,"EpipolarTiles : crop outside the source");
    }

    // Coverage of the region by the tiles, and exact equality of the assembly with the global resampling
    for (const auto & aPix : aPainted.DIm())
        MMVII_INTERNAL_ASSERT_bench(aPainted.DIm().GetV(aPix)==1, "EpipolarTiles : region not covered by the tiles");
    CmpWindow(aCanvas,aGlobIm1,cPt2di(0,0),"EpipolarTiles : assembly of the master tiles");

    for (const auto & aName : aToRemove)
        RemoveFile(aName,SVP::Yes);

}


// ============================================================
//  BenchEpipolarNoRPC : sensors with no native Z interval (conic camera).
//  Exercises cParams::mZIntv (mandatory, else user error) and
//  GenerateSensorRPC's own override.
// ============================================================

static void EpipolarNoRPCBody(cParamExeBench & aParam)
{

    const std::string aTmpDir = cMMVII_Appli::CurrentAppli().TmpDirTestMMVII() + "EpipolarNoRPC/";
    CreateDirectories(aTmpDir);

    // Synthetic conic cameras: no RPC, hence no native Z interval.
    std::unique_ptr<cPerspCamIntrCalib> aCalib(cPerspCamIntrCalib::SimpleCalib("SimulConic",cPt2di(4000,3000),4000.0));
    const cPt3dr aTarget(0.0,0.0,0.0);
    std::unique_ptr<cSensorCamPC> aCam1(BuildLookAtConicCam("Conic1",cPt3dr(-100.0,0.0,500.0),aTarget,aCalib.get()));
    std::unique_ptr<cSensorCamPC> aCam2(BuildLookAtConicCam("Conic2",cPt3dr( 100.0,0.0,500.0),aTarget,aCalib.get()));

    MMVII_INTERNAL_ASSERT_bench(! aCam1->HasIntervalZ(), "Synthetic conic camera unexpectedly has a Z interval");
    MMVII_INTERNAL_ASSERT_bench(! aCam2->HasIntervalZ(), "Synthetic conic camera unexpectedly has a Z interval");

    // Ground altitude range around the target plane (z=0)
    const cPt2dr aZIntv(-50.0,50.0);

    // ---- Missing ZIntv, no native interval : must raise a user error
    {
        bool aGotExpectedError = false;
        {
            cBenchErrorCatcher aCatcher;
            try
            {
                auto aParams = cEpipolarRectification::cParams{3,7,3};
                cEpipolarRectification(*aCam1,*aCam2,aParams).Compute();
            }
            catch (const cBenchNoZIntvError &)
            {
                aGotExpectedError = true;
            }
        }
        MMVII_INTERNAL_ASSERT_bench(aGotExpectedError, "Expected error when ZIntv is missing and sensor has no Z interval");
    }

    // ---- Images that do not overlap : explicit user error
    {
        std::unique_ptr<cSensorCamPC> aCamFar(BuildLookAtConicCam("ConicFar",cPt3dr(49900.0,0.0,500.0),cPt3dr(50000.0,0.0,0.0),aCalib.get()));
        bool aGotExpectedError = false;
        {
            cBenchErrorCatcher aCatcher;
            try
            {
                auto aParamsFar = cEpipolarRectification::cParams{3,7,3};
                aParamsFar.mZIntv = aZIntv;
                aParamsFar.mNoWarnings = true;
                cEpipolarRectification(*aCam1,*aCamFar,aParamsFar).Compute();
            }
            catch (const cBenchNoZIntvError &)
            {
                aGotExpectedError = true;
            }
        }
        MMVII_INTERNAL_ASSERT_bench(aGotExpectedError, "Expected error for two images that do not overlap");
        MMVII_INTERNAL_ASSERT_bench(TheBenchLastErrorMes.find("do not overlap") != std::string::npos, "Overlap error : unexpected message " + TheBenchLastErrorMes);
    }

    // ---- ZIntv provided : rectification succeeds
    auto aParams = cEpipolarRectification::cParams{3,7,3};
    aParams.mZIntv = aZIntv;
    auto aRectifier = cEpipolarRectification(*aCam1,*aCam2,aParams);
    auto aEpipModel = aRectifier.Compute();
    MMVII_INTERNAL_ASSERT_bench(aRectifier.NbPairs12() > 0, "No H-compatible pairs with ZIntv override");
    MMVII_INTERNAL_ASSERT_bench(aRectifier.NbPairs21() > 0, "No H-compatible pairs with ZIntv override");

    // Also exercise GenerateSensorRPC with the same override; files removed at the end.
    auto aEpipRPCName1 = aTmpDir + "Epip-Conic1.xml";
    auto aResampSI1 = std::unique_ptr<cSensorImage>(aCam1->GenerateSensorRPC(&aEpipModel.EpipMap1(), nullptr, false, "Epip-Conic1", aZIntv));
    aResampSI1->ToFile(aEpipRPCName1);
    auto aEpipSensor1 = std::unique_ptr<cSensorImage>(ReadExternalSensor(aEpipRPCName1, "Epip-Conic1", false));
    MMVII_INTERNAL_ASSERT_bench(aEpipSensor1->HasIntervalZ(), "Generated epipolar RPC has no Z interval");

    RemoveRecurs(aTmpDir,true);

    // Guard of the RPC fit : abandon a grid whose max residual is too high and try the next one, then fail.
    // The threshold is calibrated on the residuals of two grids so that the test does not depend on measured values.
    {
        const auto aSavedGrids = TheRPCFitGrids;
        const tREAL8 aSavedThr = TheRPCFitMaxResPx;
        auto FitWith = [&](const std::vector<std::pair<int,int>> & aGrids,tREAL8 aThr)
        {
            TheRPCFitGrids = aGrids;
            TheRPCFitMaxResPx = aThr;
            return std::unique_ptr<cSensorImage>(aCam1->GenerateSensorRPC(&aEpipModel.EpipMap1(), nullptr, false, "Epip-Guard", aZIntv));
        };
        FitWith({{8,3}},1e9);
        const tREAL8 aRes1 = TheRPCFitLastMaxRes;
        FitWith({{16,9}},1e9);
        const tREAL8 aRes2 = TheRPCFitLastMaxRes;
        MMVII_INTERNAL_ASSERT_bench((aRes1>0) && (aRes2>0), "RPC fit guard : residuals not measured");
        const int aNbRetry0 = TheRPCFitNbRetry;
        if (aRes2 < aRes1)   // threshold between the two : the first grid is abandoned, the second one accepted
        {
            auto aSI = FitWith({{8,3},{16,9}},(aRes1+aRes2)/2);
            MMVII_INTERNAL_ASSERT_bench((TheRPCFitNbRetry==aNbRetry0+1) && aSI, "RPC fit guard : expected one retry on the next grid");
        }
        // Threshold below every residual : failure
        bool aGotExpectedError = false;
        {
            cBenchErrorCatcher aCatcher;
            try
            {
                FitWith({{8,3},{16,9}},std::min(aRes1,aRes2)/2);
            }
            catch (const cBenchNoZIntvError &)
            {
                aGotExpectedError = true;
            }
        }
        MMVII_INTERNAL_ASSERT_bench(aGotExpectedError, "RPC fit guard : expected an error when no grid is accepted");
        TheRPCFitGrids = aSavedGrids;
        TheRPCFitMaxResPx = aSavedThr;
    }

}


// ============================================================
//  BenchEpipolarCompare : closed form vs generic algorithm on the same synthetic pair of central perspective cameras.
//  The two geometries are different images, so they are compared by the property they share : the same ground point
//  has the same row in both images of the pair, through each model.
// ============================================================

static void EpipolarCompareBody(cParamExeBench & aParam)
{

    std::unique_ptr<cPerspCamIntrCalib> aCalib(cPerspCamIntrCalib::SimpleCalib("SimulConicCmp",cPt2di(4000,3000),4000.0));
    const cPt3dr aTarget(0.0,0.0,0.0);
    std::unique_ptr<cSensorCamPC> aCam1(BuildLookAtConicCam("ConicCmp1",cPt3dr(-100.0,0.0,500.0),aTarget,aCalib.get()));
    std::unique_ptr<cSensorCamPC> aCam2(BuildLookAtConicCam("ConicCmp2",cPt3dr( 100.0,0.0,500.0),aTarget,aCalib.get()));
    const cPt2dr aZIntv(-50.0,50.0);

    auto aParamsGen = cEpipolarRectification::cParams{3,7,3};
    aParamsGen.mZIntv = aZIntv;
    aParamsGen.mNoWarnings = true;
    const auto aModelGen = cEpipolarRectification(*aCam1,*aCam2,aParamsGen).Compute();

    cEpipolarRectificationPC::cParams aParamsPC;
    aParamsPC.mZIntv = aZIntv;
    aParamsPC.mNoWarnings = true;
    const auto aModelPC = cEpipolarRectificationPC(*aCam1,*aCam2,aParamsPC).Compute();

    // Same ground points (random Z in the interval) through both models
    tREAL8 aMaxRowGen = 0, aMaxRowPC = 0, aSumSqGen = 0;
    int aNbTested = 0;
    for (int aTry=0 ; (aNbTested<300) && (aTry<20000) ; aTry++)
    {
        bool isOk = false;
        const cHomogCpleIm aCple = aCam1->RandomVisibleCple(RandInInterval(aZIntv),*aCam2,10000,&isOk);
        if (! isOk)
            continue;
        aNbTested++;
        const tREAL8 aRowGen = std::abs(aModelGen.EpipMap1().Value(aCple.mP1).y() - aModelGen.EpipMap2().Value(aCple.mP2).y());
        aMaxRowGen = std::max(aMaxRowGen,aRowGen);
        aSumSqGen += Square(aRowGen);
        aMaxRowPC  = std::max(aMaxRowPC ,std::abs(aModelPC .EpipMap1().Value(aCple.mP1).y() - aModelPC .EpipMap2().Value(aCple.mP2).y()));
    }
    MMVII_INTERNAL_ASSERT_bench(aNbTested==300,"EpipolarCompare : not enough ground points");
    MMVII_INTERNAL_ASSERT_bench(aMaxRowPC<1e-3,"EpipolarCompare : closed form rows differ : " + ToStr(aMaxRowPC));
    // generic : fitted polynomials (low degrees here), the command bounds the sigma of the residuals by MaxResid (0.1)
    MMVII_INTERNAL_ASSERT_bench(std::sqrt(aSumSqGen/aNbTested)<0.1,"EpipolarCompare : generic rows differ (rms) : " + ToStr(std::sqrt(aSumSqGen/aNbTested)));
    MMVII_INTERNAL_ASSERT_bench(aMaxRowGen<0.5,"EpipolarCompare : generic rows differ (max) : " + ToStr(aMaxRowGen));

    // Both frames are usable and of the same order
    for (const auto & aPair : {std::make_pair(&aModelGen.EpipMap1(),&aModelPC.EpipMap1()),std::make_pair(&aModelGen.EpipMap2(),&aModelPC.EpipMap2())})
    {
        const cPt2di aSzGen = aPair.first->EpipImSz(), aSzPC = aPair.second->EpipImSz();
        MMVII_INTERNAL_ASSERT_bench((aSzGen.x()>0) && (aSzPC.x()>0) && (aSzGen.y()>0) && (aSzPC.y()>0),"EpipolarCompare : empty frame");
        MMVII_INTERNAL_ASSERT_bench((aSzPC.x()<2*aSzGen.x()) && (aSzGen.x()<2*aSzPC.x()) && (aSzPC.y()<2*aSzGen.y()) && (aSzGen.y()<2*aSzPC.y()),
                                    "EpipolarCompare : frame sizes of the two algorithms differ by more than a factor 2");
    }

}

// ============================================================
//  BenchEpipolarZFromTieP : Z inferred from tie points, and priority order
//  (ZIntv must win over a valid TieP-derived interval).
// ============================================================

static void EpipolarZFromTiePBody(cParamExeBench & aParam)
{

    std::unique_ptr<cPerspCamIntrCalib> aCalib(cPerspCamIntrCalib::SimpleCalib("SimulConicZT",cPt2di(4000,3000),4000.0));
    const cPt3dr aTarget(0.0,0.0,0.0);
    std::unique_ptr<cSensorCamPC> aCam1(BuildLookAtConicCam("ConicZT1",cPt3dr(-100.0,0.0,500.0),aTarget,aCalib.get()));
    std::unique_ptr<cSensorCamPC> aCam2(BuildLookAtConicCam("ConicZT2",cPt3dr( 100.0,0.0,500.0),aTarget,aCalib.get()));

    // Synthetic tie points at a random Z within a known range, standing in for
    // real matched points.
    const cPt2dr aKnownZ(-30.0,30.0);
    cSetHomogCpleIm aSetH;
    int aTries = 0;
    while ((aSetH.NbH() < 300) && (aTries < 20000))
    {
        ++aTries;
        bool isOk = false;
        cHomogCpleIm aCple = aCam1->RandomVisibleCple(RandInInterval(aKnownZ),*aCam2,10000,&isOk);
        if (isOk)
            aSetH.Add(aCple);
    }
    MMVII_INTERNAL_ASSERT_bench(aSetH.NbH() >= 300, "Could not generate enough synthetic tie points");

    // ---- TieP alone : Z interval inferred close to the known range, succeeds
    {
        auto aParams = cEpipolarRectification::cParams{3,7,3};
        aParams.mHomolPts = aSetH;
        auto aRectifier = cEpipolarRectification(*aCam1,*aCam2,aParams);
        auto anEpipModel = aRectifier.Compute();
        MMVII_INTERNAL_ASSERT_bench(aRectifier.NbPairs12() > 0, "No H-compatible pairs with TieP-derived Z");
        MMVII_INTERNAL_ASSERT_bench(aRectifier.NbPairs21() > 0, "No H-compatible pairs with TieP-derived Z");

        const cPt2dr aUsed = anEpipModel.EpipMap1().ZInterval();
        MMVII_INTERNAL_ASSERT_bench((aUsed.x() <= aKnownZ.x()) && (aUsed.x() > aKnownZ.x()-20.0),
            "TieP-derived Zmin implausible : " + ToStr(aUsed.x()));
        MMVII_INTERNAL_ASSERT_bench((aUsed.y() >= aKnownZ.y()) && (aUsed.y() < aKnownZ.y()+20.0),
            "TieP-derived Zmax implausible : " + ToStr(aUsed.y()));
    }

    // ---- ZIntv + TieP both given : ZIntv must win (and warn). Checked directly via
    // ZIntervalUsed1/2, since Compute() could succeed either way here.
    {
        const cPt2dr anAbsurdZIntv(1.0e6,1.0e6 + 1.0);
        auto aParams = cEpipolarRectification::cParams{3,7,3};
        aParams.mHomolPts = aSetH;
        aParams.mZIntv = anAbsurdZIntv;
        aParams.mNoWarnings = true;  // suppress the expected warning about ZIntv overriding TieP-derived Z
        auto aRectifier = cEpipolarRectification(*aCam1,*aCam2,aParams);
        auto aEpipModel = aRectifier.Compute();
        MMVII_INTERNAL_ASSERT_bench(aEpipModel.EpipMap1().ZInterval() == anAbsurdZIntv,
            "ZIntv did not take priority over TieP-derived Z for camera 1");
        MMVII_INTERNAL_ASSERT_bench(aEpipModel.EpipMap2().ZInterval() == anAbsurdZIntv,
            "ZIntv did not take priority over TieP-derived Z for camera 2");
    }

    // ---- Closed form refused for a 360 degrees (EquiRect) camera, with an explicit error
    {
        std::unique_ptr<cPerspCamIntrCalib> aCalibEq(cPerspCamIntrCalib::RandomCalib(eProjPC::eEquiRect,0));
        cSensorCamPC aEq1("EquiRect1",cIsometry3D<tREAL8>::Identity(),aCalibEq.get());
        cSensorCamPC aEq2("EquiRect2",cIsometry3D<tREAL8>(cPt3dr(1,0,0),cRotation3D<tREAL8>::Identity()),aCalibEq.get());
        bool aGotExpectedError = false;
        {
            cBenchErrorCatcher aCatcher;
            try
            {
                cEpipolarRectificationPC::cParams aParams;
                aParams.mZIntv = cPt2dr(0,1);
                cEpipolarRectificationPC(aEq1,aEq2,aParams).Compute();
            }
            catch (const cBenchNoZIntvError &)
            {
                aGotExpectedError = true;
            }
        }
        MMVII_INTERNAL_ASSERT_bench(aGotExpectedError && (TheBenchLastErrorMes.find("EquiRect") != std::string::npos),"Closed form : expected an explicit error for EquiRect cameras");
    }

    // ---- Closed form (central perspective cameras) : same rules for the Z interval
    {
        cEpipolarRectificationPC::cParams aParams;
        aParams.mHomolPts = aSetH;
        const auto aModelTieP = cEpipolarRectificationPC(*aCam1,*aCam2,aParams).Compute();
        for (const auto & aMap : {&aModelTieP.EpipMap1(),&aModelTieP.EpipMap2()})
        {
            const cPt2dr aUsed = aMap->ZInterval();
            MMVII_INTERNAL_ASSERT_bench((aUsed.x() <= aKnownZ.x()) && (aUsed.x() > aKnownZ.x()-20.0),"Closed form : TieP-derived Zmin implausible : " + ToStr(aUsed.x()));
            MMVII_INTERNAL_ASSERT_bench((aUsed.y() >= aKnownZ.y()) && (aUsed.y() < aKnownZ.y()+20.0),"Closed form : TieP-derived Zmax implausible : " + ToStr(aUsed.y()));
        }

        const cPt2dr anAbsurdZIntv(1.0e6,1.0e6 + 1.0);
        aParams.mZIntv = anAbsurdZIntv;
        aParams.mNoWarnings = true;
        const auto aModelZIntv = cEpipolarRectificationPC(*aCam1,*aCam2,aParams).Compute();
        MMVII_INTERNAL_ASSERT_bench((aModelZIntv.EpipMap1().ZInterval() == anAbsurdZIntv) && (aModelZIntv.EpipMap2().ZInterval() == anAbsurdZIntv),"Closed form : ZIntv did not take priority over TieP-derived Z");
    }

    // ---- Too few tie points after filtering : error, not a silent small interval.
    {
        bool aGotExpectedError = false;
        {
            cBenchErrorCatcher aCatcher;
            try
            {
                auto aParams = cEpipolarRectification::cParams{3,7,3};
                aParams.mHomolPts = aSetH;
                aParams.mTiePMaxRes = -1.0;
                cEpipolarRectification(*aCam1,*aCam2,aParams).Compute();
            }
            catch (const cBenchNoZIntvError &)
            {
                aGotExpectedError = true;
            }
        }
        MMVII_INTERNAL_ASSERT_bench(aGotExpectedError,
            "Expected error when too few tie points survive residual filtering");
    }

}


// ============================================================
//  Groups of epipolar benches (one registered bench per topic)
// ============================================================

void BenchEpipolarCrop(cParamExeBench & aParam)
{
    if (! aParam.NewBench("EpipolarCrop")) return;
    EpipolarResamplingBody(aParam);   // crop and validity mask of the resampling
    EpipolarSlaveCropBody(aParam);    // slave crop derived from a master crop and the Z interval
    EpipolarTilesBody(aParam);        // tiles of a pair equal the global resampling
    aParam.EndBench();
}

void BenchEpipolarZ(cParamExeBench & aParam)
{
    if (! aParam.NewBench("EpipolarZ")) return;
    EpipolarNoRPCBody(aParam);        // sensors with no native Z interval
    EpipolarZFromTiePBody(aParam);    // Z interval from tie points, priority of ZIntv, closed form, EquiRect refused
    aParam.EndBench();
}

void BenchEpipolarPC(cParamExeBench & aParam)
{
    if (! aParam.NewBench("EpipolarPC")) return;
    BenchEpipolarPCBody(aParam);      // closed form: rows, round trips, virtual cameras, model on disk
    EpipolarCompareBody(aParam);      // closed form vs generic on the same pair
    aParam.EndBench();
}

} // namespace MMVII
