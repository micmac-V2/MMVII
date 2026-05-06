#include "cEpipolarRectification.h"
#include "MMVII_Error.h"
#include "MMVII_Geom2D.h"
#include "MMVII_Sensor.h"
#include <cmath>
#include <cassert>

namespace MMVII {


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
    return cPt2dr(q.x(), mV.Eval(q));
}

cPt2dr cEpipPolyMapping::Inverse(const cPt2dr& aPt) const
{
    return FromRotatedFrame(cPt2dr(aPt.x(), mW.Eval(aPt)));
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

cEpipPolyModel cEpipolarRectification::Compute() const
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
    aDir2 = - aDir2;

    // TODOCM: Check aDir1 ~= 0 0 && aDir2 ~= 0
    // if ((dir2.x+dir1.x) <0)
    // {
    //     dir1 = -dir1;
    //     dir2 = -dir2;
    // }
    // TODOCM: Check d1 and d2 /= 0 ?
    aDir1 = VUnit(aDir1);
    aDir2 = VUnit(aDir2);

    StdOut() << "Nb Pairs 1->2 : " << aPairsA.size() << std::endl;
    StdOut() << "Nb Pairs 2->1 : " << aPairsB.size() << std::endl;
    StdOut() << "C1: " << aCenter1 << " Dir1:" << aDir1 << std::endl;
    StdOut() << "C2: " << aCenter2 << " Dir2:" << aDir2 << std::endl;

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

    return cEpipPolyModel {
        std::make_unique<cEpipPolyMapping>(aV1,aW1,aCenter1,aDir1),
        std::make_unique<cEpipPolyMapping>(aV2,aW2,aCenter2,aDir2),
    };
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
        cPolyXY_Nd&           aV2) const
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
    StdOut() << "V1,V2 var = " << aSolver.VarCurSol() << std::endl;

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
        UseFromPair                  aUsePt) const
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
    StdOut() << "W" << (aUsePt == UseFromPair::PT1 ? '1' : '2') << " var = " <<  aSolver.VarCurSol() << std::endl;

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
    const double aEpsMarginRel = mParams.mEpsMarginRel;
    // Lambda to convert Z steps to "world" Z
    auto Step2Z = [&](int aKZ) -> tREAL8 {
        tREAL8 aW = (aKZ + aEpsMarginRel) / (nZ - 1 + 2 * aEpsMarginRel);
        return (Zmin * (1 - aW)) + Zmax * aW;
    };

    std::vector<cPt2dr> aVPts = aCamM.PtsSampledOnSensor(nXY, aEpsMarginRel);

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
                aCamM.ImageAndZ2Ground(cPt3dr(pM.x(), pM.y(), Z0));
            const cPt3dr aGround1 =
                aCamM.ImageAndZ2Ground(cPt3dr(pM.x(), pM.y(), Z1));

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
    const cDataGenUnTypedIm<2> *aIm1,
    const cDataGenUnTypedIm<2> *aIm2,
    eEpipFrm aFrmType,
    int aMargin
    )
{

    // TODOCM : use margin
    // TODOCM: add option : minimum, maximum, im1 , im2
    auto frame1 = EpipMap1().BoxOfFrontier(aIm1->ToR(),1.0);
    auto frame2 = EpipMap2().BoxOfFrontier(aIm2->ToR(),1.0);

    auto P1_0 = frame1.P0()-cPt2dr(aMargin,aMargin);
    auto P1_1 = frame1.P1()+cPt2dr(aMargin,aMargin);;
    auto P2_0 = frame2.P0()-cPt2dr(aMargin,aMargin);;
    auto P2_1 = frame2.P1()+cPt2dr(aMargin,aMargin);;;

    double yMin=0;
    double yMax=0;
    FakeUseIt(yMin);
    FakeUseIt(yMax);
    switch(aFrmType) {
    case eEpipFrm::eNbVals:
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
    }

    P1_0.y() = P2_0.y() = yMin;
    P1_1.y() = P2_1.y() = yMax;

    mFrame1 = cTplBox<double,2>(P1_0,P1_1).ToI();
    mFrame2 = cTplBox<double,2>(P2_0,P2_1).ToI();
}

} // namespace MMVII
