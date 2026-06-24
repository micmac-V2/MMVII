#include "cEpipolarRectification.h"
#include "MMVII_Error.h"
#include "MMVII_Geom2D.h"
#include "MMVII_Sensor.h"
#include "../Sensors/cExternalSensor.h"

#include <cmath>
#include <cassert>

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

    // Altitude range from the master camera's RPC validity interval
    const cPt2dr aZInterval = aCamM.GetIntervalZ();
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
    for (const auto& [aES1,aES2] : {std::make_pair(&aEpipSensor1, &aEpipSensor2),std::make_pair(&aEpipSensor2, &aEpipSensor1)})
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


} // namespace MMVII
