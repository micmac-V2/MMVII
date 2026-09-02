#include "cEpipolarRectification.h"
#include "MMVII_Error.h"
#include "MMVII_Geom2D.h"
#include "MMVII_Sensor.h"
#include "MMVII_PCSens.h"
#include "../Sensors/cExternalSensor.h"

#include <cmath>
#include <cassert>
#include <algorithm>

namespace MMVII {

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

    std::vector<cEpiPair> aPairsA, aPairsB;
    cPt2dr aCenter1, aCenter2;
    cPt2dr aDir1,    aDir2;

    GenerateData(mCam1, mCam2, aPairsA, aCenter1, aDir2);
    GenerateData(mCam2, mCam1, aPairsB, aCenter2, aDir1);

    // We must inverse aDir2 because it is computed in the direction from I2 to I1, but we want it in the direction from I1 to I2
    aDir2 = - aDir2;

    if ((aDir2.x() + aDir1.x()) <0)
    {
        aDir1 = -aDir1;
        aDir2 = -aDir2;
    }

    // TODOCM: Check d1 and d2 /= 0 ?
    aDir1 = VUnit(aDir1);
    aDir2 = VUnit(aDir2);
    mNbPairs12 = aPairsA.size();
    mNbPairs21 = aPairsB.size();

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

    // aPairsA already stores (masterPt in I1, slavePt in I2)
    // aPairsB stores          (masterPt in I2, slavePt in I1)
    //   => swap to keep the convention (pt1, pt2)

    std::vector<cEpiPair> aRotPairs;
    aRotPairs.reserve(aPairsA.size() + aPairsB.size());

    for (const auto& pr : aPairsA)
        aRotPairs.push_back({ Rotate(pr.mP1, aCenter1, aDir1),
                              Rotate(pr.mP2, aCenter2, aDir2) });

    for (const auto& pr : aPairsB)
        aRotPairs.push_back({ Rotate(pr.mP2, aCenter1, aDir1),   // I1 pt
                              Rotate(pr.mP1, aCenter2, aDir2) }); // I2 pt

    // ----------------------------------------------------------
    //  Step 3 – estimate V1 (with Y-axis identity) and V2
    // ----------------------------------------------------------
    cPolyXY_N<tREAL8> aV1(mParams.mPolyDegree);
    cPolyXY_N<tREAL8> aV2(mParams.mPolyDegree);
    EstimateForwardPolynomials(aRotPairs, aV1, aV2);

    // ----------------------------------------------------------
    //  Step 4 – estimate inverse polynomials W1, W2
    // ----------------------------------------------------------

    cPolyXY_N<tREAL8> aW1(mParams.mPolyDegreeInv);
    cPolyXY_N<tREAL8> aW2(mParams.mPolyDegreeInv);
    EstimateInversePolynomial(aRotPairs, aV1, aW1, UseFromPair::PT1);
    EstimateInversePolynomial(aRotPairs, aV2, aW2, UseFromPair::PT2);

    auto anEpipPolyModel = cEpipPolyModel {
        std::make_unique<cEpipPolyMapping>(aV1,aW1,aCenter1,aDir1),
        std::make_unique<cEpipPolyMapping>(aV2,aW2,aCenter2,aDir2),
    };
    anEpipPolyModel.ComputeCommonFraming(mCam1.PixelDomain().Box(),mCam2.PixelDomain().Box(),mParams.mEpipFrm, mParams.mMargin);

    // TODOCM: Prendre 1 pt sur 2 pour calcul et les autres pour la verif

    return anEpipPolyModel;
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

        // RHS = -locked contribution of V1 at q1 = -q1.y()
        const double aRHS = -aV1IdOnY.LockedContribution(q1);

        aSolver.PublicAddObservation(1.0, aCoeff, aRHS);
    }

    const cDenseVect<double> aSol = aSolver.PublicSolve();
    // TODOCM : error on too big variance. Scale of variance ?
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

    for (const auto& pr : aPairs)
    {
        const cPt2dr& qk = aUsePt == UseFromPair::PT1 ? pr.mP1 : pr.mP2;

        // Epipolar coordinates of qk
        const double u = qk.x();
        const double v = aVk.Eval(qk);   // v = Vk(qk)

        // Observation : Wk(u, v) = qk.y
        const cDenseVect<double> aCoeff = aWk.BasisVector(u, v);
        const double             aRHS   = qk.y();

        aSolver.PublicAddObservation(1.0, aCoeff, aRHS);
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



void cEpipolarRectification::GenerateData(const cSensorImage &aCamM,
                                          const cSensorImage &aCamS,
                                          std::vector<cEpiPair> &aOutPairs,
                                          cPt2dr &aOutCenterM,
                                          cPt2dr &aOutDirS) const {
    aOutPairs.clear();
    aOutCenterM = cPt2dr(0, 0);
    aOutDirS = cPt2dr(0, 0);

    // Steps on sensors
    const double nXY = mParams.mNbXYSteps;

    // Altitude range for this master camera : mZIntv > tie-point-derived > native.
    const cPt2dr aZInterval = EffectiveZInterval(aCamM);
    const double Zmin = aZInterval.x();
    const double Zmax = aZInterval.y();

    // Altitude step : NbZLevels levels => (NbZLevels-1) intervals
    const int nZ = mParams.mNbZLevels;

    // Margin
    // Lambda to convert Z steps to "world" Z
    auto Step2Z = [&](int aKZ) -> tREAL8 {
        tREAL8 aW = (aKZ) / (nZ - 1);
        return (Zmin * (1 - aW)) + Zmax * aW;
    };

    std::vector<cPt2dr> aVPts = aCamM.PtsSampledOnSensor(nXY, 0);

    int nAccum = 0;
    for (const auto &pM : aVPts) {
        for (int aKZ = 0; aKZ < nZ; aKZ++) {
            cPt3dr aP13D;
            const double Z0 = Step2Z(aKZ);
            const double Z1 = Step2Z(aKZ + 1);

            // 3D point on bundle at altitude Z0 and Z1
            // cSensorImage::ImageAndZ2Ground expects a cPt3dr
            // with x=col, y=row, z=altitude
            const cPt3dr aGround0 =
                aCamM.ImageAndZ2Ground(TP3z(pM, Z0));
            const cPt3dr aGround1 =
                aCamM.ImageAndZ2Ground(TP3z(pM, Z1));

            // Project into the slave image
            const cPt2dr pS0 = aCamS.Ground2Image(aGround0);
            const cPt2dr pS1 = aCamS.Ground2Image(aGround1);

            if (!aCamS.IsVisibleOnImFrame(pS0))
                continue;
            if (!aCamS.IsVisibleOnImFrame(pS1))
                continue;

            // Store the pair at Z0
            aOutPairs.push_back({pM, pS0});

            // Accumulate centroid
            aOutCenterM = aOutCenterM + pM;

            // Accumulate normalised epipolar direction in slave
            cPt2dr aDelta = pS1 - pS0;
            if (SqN2(aDelta) > 1e-16) {
                aDelta = VUnit(aDelta);
            }
            aOutDirS = aOutDirS + aDelta;

            ++nAccum;
        }
    }

    MMVII_INTERNAL_ASSERT_User(
        nAccum > 50, eTyUEr::eUnClassedError,
        "Insufficient number of common points on two images");
    aOutCenterM = aOutCenterM * (1.0 / nAccum);
    aOutDirS = VUnit(aOutDirS);
}

// ============================================================
//  EffectiveZInterval : mZIntv > tie-point-derived > aCamM's own native.
//  Overriding a lower-priority source only warns.
// ============================================================

cPt2dr cEpipolarRectification::EffectiveZInterval(const cSensorImage & aCamM) const
{
    cPt2dr aResult(0,0);

    if (mParams.mZIntv)
    {
        if (aCamM.HasIntervalZ() && !mParams.mNoWarnings)
        {
            MMVII_USER_WARNING("Provided ZIntv overrides sensor's own Z validity interval");
        }
        if (mParams.mHomolPts && !mParams.mNoWarnings)
        {
            MMVII_USER_WARNING("Provided ZIntv overrides tie-point-derived Z validity interval");
        }
        aResult = *mParams.mZIntv;
    }
    else if (mParams.mHomolPts)
    {
        aResult = ZIntervalFromHomolPts();
        if (aCamM.HasIntervalZ() && !mParams.mNoWarnings)
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

    ((&aCamM == &mCam1) ? mZIntervalUsed1 : mZIntervalUsed2) = aResult;
    return aResult;
}

// ============================================================
//  ZIntervalFromHomolPts : Z envelope of triangulated tie points, inflated
//  by mZMargin. Memoized : identical for both master cameras.
// ============================================================

cPt2dr cEpipolarRectification::ZIntervalFromHomolPts() const
{
    if (mCachedHomolZIntv)
        return *mCachedHomolZIntv;

    int aNbKept = 0;
    tREAL8 aZmin = 0.0;
    tREAL8 aZmax = 0.0;
    for (const auto & aCple : mParams.mHomolPts->SetH())
    {
        const tREAL8 aRes = mCam1.PixResInterBundle(aCple, mCam2);
        if (aRes > mParams.mTiePMaxRes)
            continue;

        const tREAL8 aZ = mCam1.PInterBundle(aCple, mCam2).z();
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

    const cPt2dr aSz = mCam1.PixelDomain().Box().Sz();
    const int aMinNb = std::max(mParams.mTiePMinNbFloor,
                                 (int)std::ceil(mParams.mTiePMinNbRatio * std::sqrt(aSz.x() * aSz.y())));
    MMVII_INTERNAL_ASSERT_User(aNbKept >= aMinNb, eTyUEr::eUnClassedError,
        "Not enough tie points after residual filtering to infer Z interval (" + ToStr(aNbKept)
        + " < " + ToStr(aMinNb) + "); provide ZIntv=[Zmin,Zmax], relax TiePMaxRes/TiePMinNb*, or add tie points");

    MMVII_INTERNAL_ASSERT_User((aZmax - aZmin) > 1e-6, eTyUEr::eUnClassedError,
        "Tie points give a degenerate (near-flat) Z interval [" + ToStr(aZmin) + "," + ToStr(aZmax)
        + "]; provide ZIntv=[Zmin,Zmax] explicitly for this scene");

    const tREAL8 aMargin = mParams.mZMargin * (aZmax - aZmin);
    mCachedHomolZIntv = cPt2dr(aZmin - aMargin, aZmax + aMargin);
    return *mCachedHomolZIntv;
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
    // TODOCM : Randomize parameters (ou faire une loop sur degree)
    auto aParams = cEpipolarRectification::cParams{5,9,100,3};
    auto aRectifier = cEpipolarRectification(*aSensor1, *aSensor2, aParams);
    auto aEpipModel = aRectifier.Compute();

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


namespace {

// Catches MMVII_UserError/ASSERT_User via a throwing handler instead of abort();
// RAII-restored, same pattern as cProfileErrorCatcher.
struct cBenchNoZIntvError {};

void BenchNoZIntvErrorHandler(const std::string &, const std::string &, const char *, int)
{
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
//  BenchEpipolarNoRPC : sensors with no native Z interval (conic camera).
//  Exercises cParams::mZIntv (mandatory, else user error) and
//  GenerateSensorRPC's own override.
// ============================================================

void BenchEpipolarNoRPC(cParamExeBench & aParam)
{
    if (! aParam.NewBench("EpipolarNoRPC")) return;

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
                auto aParams = cEpipolarRectification::cParams{3,7,20,3};
                cEpipolarRectification(*aCam1,*aCam2,aParams).Compute();
            }
            catch (const cBenchNoZIntvError &)
            {
                aGotExpectedError = true;
            }
        }
        MMVII_INTERNAL_ASSERT_bench(aGotExpectedError, "Expected error when ZIntv is missing and sensor has no Z interval");
    }

    // ---- ZIntv provided : rectification succeeds
    auto aParams = cEpipolarRectification::cParams{3,7,20,3};
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

    aParam.EndBench();
}


// ============================================================
//  BenchEpipolarZFromTieP : Z inferred from tie points, and priority order
//  (ZIntv must win over a valid TieP-derived interval).
// ============================================================

void BenchEpipolarZFromTieP(cParamExeBench & aParam)
{
    if (! aParam.NewBench("EpipolarZFromTieP")) return;

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
        auto aParams = cEpipolarRectification::cParams{3,7,20,3};
        aParams.mHomolPts = aSetH;
        auto aRectifier = cEpipolarRectification(*aCam1,*aCam2,aParams);
        aRectifier.Compute();
        MMVII_INTERNAL_ASSERT_bench(aRectifier.NbPairs12() > 0, "No H-compatible pairs with TieP-derived Z");
        MMVII_INTERNAL_ASSERT_bench(aRectifier.NbPairs21() > 0, "No H-compatible pairs with TieP-derived Z");

        const cPt2dr aUsed = aRectifier.ZIntervalUsed1();
        MMVII_INTERNAL_ASSERT_bench((aUsed.x() <= aKnownZ.x()) && (aUsed.x() > aKnownZ.x()-20.0),
            "TieP-derived Zmin implausible : " + ToStr(aUsed.x()));
        MMVII_INTERNAL_ASSERT_bench((aUsed.y() >= aKnownZ.y()) && (aUsed.y() < aKnownZ.y()+20.0),
            "TieP-derived Zmax implausible : " + ToStr(aUsed.y()));
    }

    // ---- ZIntv + TieP both given : ZIntv must win (and warn). Checked directly via
    // ZIntervalUsed1/2, since Compute() could succeed either way here.
    {
        const cPt2dr anAbsurdZIntv(1.0e6,1.0e6 + 1.0);
        auto aParams = cEpipolarRectification::cParams{3,7,20,3};
        aParams.mHomolPts = aSetH;
        aParams.mZIntv = anAbsurdZIntv;
        aParams.mNoWarnings = true;  // suppress the expected warning about ZIntv overriding TieP-derived Z
        auto aRectifier = cEpipolarRectification(*aCam1,*aCam2,aParams);
        aRectifier.Compute();
        MMVII_INTERNAL_ASSERT_bench(aRectifier.ZIntervalUsed1() == anAbsurdZIntv,
            "ZIntv did not take priority over TieP-derived Z for camera 1");
        MMVII_INTERNAL_ASSERT_bench(aRectifier.ZIntervalUsed2() == anAbsurdZIntv,
            "ZIntv did not take priority over TieP-derived Z for camera 2");
    }

    // ---- Too few tie points after filtering : error, not a silent small interval.
    {
        bool aGotExpectedError = false;
        {
            cBenchErrorCatcher aCatcher;
            try
            {
                auto aParams = cEpipolarRectification::cParams{3,7,20,3};
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

    aParam.EndBench();
}


} // namespace MMVII
