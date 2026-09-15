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
{ //GetCol(mVGa,aRSEV.EigenVectors(),0);
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
}
/*
 *     _____
 * P1 *     /|
 *   /     / |
 *  /___P2*  |
 *  | P4x |  *P3   x - pt in the middle
 *  |     | /      * - pts in the corners
 *P5*_____|/
 *
 *
 * */

void cGenGauss3D::GetDistrib5Pts(std::vector<cPt3dr> & aVPts,double aScale)
{
    // Correction factor: the 5-point scheme reconstructs variance as (4/5)*v^2*lambda,
    // so scale by sqrt(5)/(2*v) to recover exact eigenvalues.
    const double v = InvErrFonc(2.0/3.0);
    double aCorr = std::sqrt(5.0)/(2.0*v);

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
    aP1 = mCDG + aFact0 * InvErrFoncRationel(-1*2,2+1)
               + aFact1 * InvErrFoncRationel(1*2,2+1)
               + aFact2 * InvErrFoncRationel(1*2,2+1);
    aVPts.push_back(aP1);

    cPt3dr aP2;
    aP2 = mCDG + aFact0 * InvErrFoncRationel(1*2,2+1)
               + aFact1 * InvErrFoncRationel(-1*2,2+1)
               + aFact2 * InvErrFoncRationel(1*2,2+1);
    aVPts.push_back(aP2);

    cPt3dr aP3;
    aP3 = mCDG + aFact0 * InvErrFoncRationel(1*2,2+1)
               + aFact1 * InvErrFoncRationel(1*2,2+1)
               + aFact2 * InvErrFoncRationel(-1*2,2+1);
    aVPts.push_back(aP3);

    cPt3dr aP4;
    aP4 = mCDG + aFact0 * InvErrFoncRationel(0,2+1)
               + aFact1 * InvErrFoncRationel(0,2+1)
               + aFact2 * InvErrFoncRationel(0,2+1);
    aVPts.push_back(aP4);

    cPt3dr aP5;
    aP5 = mCDG + aFact0 * InvErrFoncRationel(-1*2,2+1)
               + aFact1 * InvErrFoncRationel(-1*2,2+1)
               + aFact2 * InvErrFoncRationel(-1*2,2+1);
    aVPts.push_back(aP5);

}

void cGenGauss3D::GetDistrib9Pts(std::vector<cPt3dr> & aVPts,double aScale)
{

}

void cGenGauss3D::GetDistrib27Pts(std::vector<cPt3dr> & aVPts,double aScale)
{

}

void cGenGauss3D::Bench()
{
    int aNbIter=5;

    for (int aKI=0; aKI<aNbIter; aKI++)
    {
        //StdOut() << "==== Iter " << aKI << std::endl;
        int aNbPts = 5 + RandUnif_N(20);

        // weighted covariance of the 3D points
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
                MMVII_INTERNAL_ASSERT_bench(aDistVec<1e-10,"cGenGauss3D::Bench");
            }
        }
    }
}

};
