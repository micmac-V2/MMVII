#include "MMVII_Geom3D.h"

namespace MMVII
{

// value that corrects for the fact that the function is an approximation
static const double A = 0.147;

double InvErrFonc(double x)
{
    static double DePiSA = 2/(M_PI*A);

    double X2 = x * x;
    double Log1MX2 = log(1-X2);

    double aV1 = Square(DePiSA + Log1MX2/2)  - Log1MX2/A;
    aV1 = sqrt(std::max(0.0,aV1)) - (DePiSA+Log1MX2/2);

    aV1 = sqrt(std::max(0.0,aV1));


    return (x>0) ? aV1 : (-aV1);
}

double InvErrFoncRationel(int P,int Q)
{
    int aSign = (P>0 ? 1 : -1) * (Q>0 ? 1 : -1);
    P = abs(P);
    Q = abs(Q);
    static std::vector<std::vector<double>  > mBuf;

    for (int q=mBuf.size() ; q<=Q ; q++)
    {
        mBuf.push_back(std::vector<double>());
        for (int p=0 ; p<q ; p++)
            mBuf.back().push_back(InvErrFonc(p/double(q)));
    }

    return mBuf[Q][P] * aSign;
}
static double FactCorrectif(int aNb)
{
    return sqrt(2.) / (1 - 0.3333/(aNb+0.5));
}

/* *********************************************************** */
/*                                                             */
/*                         cGenGaus3D                          */
/*                                                             */
/* *********************************************************** */

cGenGauss3D::cGenGauss3D(const cDenseMatrix<double> & aVecEig,
                         const cDenseVect<double> & aValEig,
                         const cDenseVect<double> & aCG) :
    mCDG(cPt3dr(aCG(0), aCG(1), aCG(2))),
    mVP(aValEig),
    mVecP(aVecEig)
{
}

void cGenGauss3D::GetDistribGaus(std::vector<cPt3dr> & aVPts,int aN1,int aN2,int aN3)
{
    cPt3dr aVec0,aVec1,aVec2;
    GetCol(aVec0,mVecP,0);
    GetCol(aVec1,mVecP,1);
    GetCol(aVec2,mVecP,2);

    aVPts.clear();

    cPt3dr aFact1 = aVec0 * (FactCorrectif(aN1) * Sqrt(mVP(0)));
    cPt3dr aFact2 = aVec1 * (FactCorrectif(aN2) * Sqrt(mVP(1)));
    cPt3dr aFact3 = aVec2 * (FactCorrectif(aN3) * Sqrt(mVP(2)));

    for (int aK1 =-aN1 ; aK1<=aN1 ; aK1++)
    {
        for (int aK2 =-aN2 ; aK2<=aN2 ; aK2++)
        {
            for (int aK3 =-aN3 ; aK3<=aN3 ; aK3++)
            {
                cPt3dr aP  =   mCDG
                           + aFact1 * InvErrFoncRationel(2*aK1,2*aN1+1)
                           + aFact2 * InvErrFoncRationel(2*aK2,2*aN2+1)
                           + aFact3 * InvErrFoncRationel(2*aK3,2*aN3+1)
                    ;
                aVPts.push_back(aP);
            }
        }
    }

}

void cGenGauss3D::GetDistribNPts(std::vector<cPt3dr> & aVPts,eDistrVirTPs aType,double aSca)
{
    switch (aType)
    {
        case eDistrVirTPs::e5Pts : GetDistrib5Pts(aVPts,aSca); break;
        case eDistrVirTPs::e9Pts : GetDistrib9Pts(aVPts,aSca); break;
        case eDistrVirTPs::e27Pts : GetDistrib27Pts(aVPts,aSca); break;
        default : MMVII_INTERNAL_ERROR("cGenGauss3D::GetDistribNPts: unhandled eDistrVirTPs value");
    }

    MMVII_INTERNAL_ASSERT_medium((int)aVPts.size()==NbPtsOfDistrib(aType),"cGenGauss3D::GetDistribNPts: wrong number of points");
}

/*
 *     _____
 * P1 *     /|
 *   /     / |
 *  /___P2*  |
 *  | P4x |  *P3   x - center of gravity point
 *  |     | /      * - the remaining points
 *P5*_____|/
 *
 *
 * */

void cGenGauss3D::GetDistrib5Pts(std::vector<cPt3dr> & aVPts,double aScale)
{
    // 4 corners at (−,+,+), (+,−,+), (+,+,−), centre, (−,−,−) : per-axis variance is
    // (4/5)*a^2*lambda, so a = sqrt(5)/2 reproduces the eigenvalues exactly
    double aCorr = std::sqrt(5.0)/(2.0);

    if (aScale<1.0)
        aCorr *= aScale;

    cPt3dr aVec0,aVec1,aVec2;
    GetCol(aVec0,mVecP,0);
    GetCol(aVec1,mVecP,1);
    GetCol(aVec2,mVecP,2);

    cPt3dr aFact0 = aVec0 * (aCorr * Sqrt(mVP(0)));
    cPt3dr aFact1 = aVec1 * (aCorr * Sqrt(mVP(1)));
    cPt3dr aFact2 = aVec2 * (aCorr * Sqrt(mVP(2)));

    aVPts.clear();

    cPt3dr aP1;
    aP1 = mCDG + aFact0 * (-1.0)
               + aFact1 * 1.0
               + aFact2 * 1.0;
    aVPts.push_back(aP1);

    cPt3dr aP2;
    aP2 = mCDG + aFact0 * 1.0
               + aFact1 * (-1.0)
               + aFact2 * 1.0;
    aVPts.push_back(aP2);

    cPt3dr aP3;
    aP3 = mCDG + aFact0 * 1.0
               + aFact1 * 1.0
               + aFact2 * (-1.0);
    aVPts.push_back(aP3);

    aVPts.push_back(mCDG);

    cPt3dr aP5;
    aP5 = mCDG + aFact0 * (-1.0)
               + aFact1 * (-1.0)
               + aFact2 * (-1.0);
    aVPts.push_back(aP5);

}

/**
 *     _____
 *    *     *|
 *   /|    / |
 *  *_|___*  |
 *  | *__x|__*     x - center of gravity (CoG) point
 *  |/    |  /      * - pts symmetric wrt (CoG)
 *  *_____* /

**/
void cGenGauss3D::GetDistrib9Pts(std::vector<cPt3dr> & aVPts,double aScale)
{
    const std::vector<int> aRange = {-1,1};

    // 8 corners at +-a*sigma on each axis plus the centre : per-axis variance is
    // (8/9)*a^2*lambda, so a = 3/sqrt(8) reproduces the eigenvalues exactly
    double aCorr = 3.0/std::sqrt(8.0);

    if (aScale<1.0)
        aCorr *= aScale;

    cPt3dr aVec0,aVec1,aVec2;
    GetCol(aVec0,mVecP,0);
    GetCol(aVec1,mVecP,1);
    GetCol(aVec2,mVecP,2);

    aVPts.clear();

    cPt3dr aFact1 = aVec0 * (aCorr * Sqrt(mVP(0)));
    cPt3dr aFact2 = aVec1 * (aCorr * Sqrt(mVP(1)));
    cPt3dr aFact3 = aVec2 * (aCorr * Sqrt(mVP(2)));

    for (auto aK1 : aRange)
    {
        for (auto aK2 : aRange)
        {
            for (auto aK3 : aRange)
            {
                cPt3dr aP  =   mCDG
                            + aFact1 * (double)aK1
                            + aFact2 * (double)aK2
                            + aFact3 * (double)aK3
                    ;
                aVPts.push_back(aP);
            }
        }
    }

    // add center of gravity
    aVPts.push_back(mCDG);
}

/**
      ______________
     *        *        *
    *        *        *
   *________*________*
     *        *        *
    *        *        *             * - generated points
   *________*________*
     *        *        *
    *        *        *
   *________*________*

**/
void cGenGauss3D::GetDistrib27Pts(std::vector<cPt3dr> & aVPts,double aScale)
{
    const int aN = 1;

    // 3x3x3 grid at {-a,0,+a}*sigma on each axis : 18 of the 27 points have a
    // non-zero coordinate on a given axis, so the per-axis variance is
    // (18/27)*a^2*lambda = (2/3)*a^2*lambda, and a = sqrt(3/2) is exact
    double aCorr = std::sqrt(1.5);

    if (aScale<1.0)
        aCorr *= aScale;

    cPt3dr aVec0,aVec1,aVec2;
    GetCol(aVec0,mVecP,0);
    GetCol(aVec1,mVecP,1);
    GetCol(aVec2,mVecP,2);

    aVPts.clear();

    cPt3dr aFact1 = aVec0 * (aCorr * Sqrt(mVP(0)));
    cPt3dr aFact2 = aVec1 * (aCorr * Sqrt(mVP(1)));
    cPt3dr aFact3 = aVec2 * (aCorr * Sqrt(mVP(2)));

    for (int aK1 =-aN ; aK1<=aN ; aK1++)
    {
        for (int aK2 =-aN ; aK2<=aN ; aK2++)
        {
            for (int aK3 =-aN ; aK3<=aN ; aK3++)
            {
                cPt3dr aP  =   mCDG
                            + aFact1 * (double)aK1
                            + aFact2 * (double)aK2
                            + aFact3 * (double)aK3
                    ;
                aVPts.push_back(aP);
            }
        }
    }
}

void cGenGauss3D::Bench()
{
    int aNbIter=5;

    for (int aKI=0; aKI<aNbIter; aKI++)
    {
        //StdOut() << "==== Iter " << aKI << std::endl;
        int aNbPts = 5 + RandUnif_N(20);

        // weighted covariance  with all weights 1.0
        cStrStat2<tREAL8> aCovMat(3);

        cPt3dr aC0 (RandUnif_C(),RandUnif_C(),RandUnif_C());
        cPt3dr aU0 (RandUnif_C(),RandUnif_C(),RandUnif_C());
        cPt3dr aU1 (RandUnif_C(),RandUnif_C(),RandUnif_C());
        cPt3dr aU2 (RandUnif_C(),RandUnif_C(),RandUnif_C());
        //StdOut() << aC0 << " " << aU0 << " " << aU1 << " " << aU2 << " " << aNbPts << std::endl;
        //getchar();
        for (int aK=0; aK<aNbPts; aK++)
        {
            cPt3dr aP = aC0 + aU0 * RandUnif_C() + aU1 * RandUnif_C() + aU2 * RandUnif_C();
            aP.x() *= 10;
            aP.y() *= 10;
            aP.z() *= 10;
            aCovMat.WeightedAdd(aP.ToVect(),1.0);
        }
        aCovMat.Normalise();

        cResulSymEigenValue<tREAL8> aVp = aCovMat.DoEigen();

        cGenGauss3D aG3D(aVp.EigenVectors(),aVp.EigenValues(),aCovMat.Moy());

        for (int aK=0; aK<int(eDistrVirTPs::eNbVals); aK++)
        {
            eDistrVirTPs aDistrib = eDistrVirTPs(aK);

            std::vector<cPt3dr> aVPts;
            aG3D.GetDistribNPts(aVPts,aDistrib,1.0);

            if (aVPts.empty()) continue; // distribution not implemented yet

            cStrStat2<tREAL8> aCovMat2(3);
            for (const auto & aP : aVPts)
                aCovMat2.WeightedAdd(aP.ToVect(),1.0);
            aCovMat2.Normalise();

            cResulSymEigenValue<tREAL8> aVp2 = aCovMat2.DoEigen();

            cGenGauss3D aG3D2(aVp2.EigenVectors(),aVp2.EigenValues(),aCovMat2.Moy());

            for(int aKV=0; aKV<3; aKV++)
            {
                double aDistVec = Sqrt(Square(aG3D.VecP(aKV)(0)-aG3D2.VecP(aKV)(0)) +
                                       Square(aG3D.VecP(aKV)(1)-aG3D2.VecP(aKV)(1)) +
                                       Square(aG3D.VecP(aKV)(2)-aG3D2.VecP(aKV)(2)));
                MMVII_INTERNAL_ASSERT_bench(aDistVec<1e-10,"cGenGauss3D::Bench VecP");

                MMVII_INTERNAL_ASSERT_bench(std::abs(aG3D2.ValP(aKV)/aG3D.ValP(aKV)-1.0)<1e-8,"cGenGauss3D::Bench ValP");

            }

            MMVII_INTERNAL_ASSERT_bench(Norm2(aG3D2.CDG()-aG3D.CDG())<1e-8,"cGenGauss3D::Bench CDG");
        }
    }
}

};
