/** \file SelecTripletsCoverage.cpp
    \brief Selection of triplets by image coverage and B/H, port of MMV1 NO_GenTripl
*/
#include <algorithm>
#include <cmath>

#include "SelecTripletsCoverage.h"
#include "MMVII_Bench.h"


namespace MMVII
{

/* ******************************** */
/*                                  */
/*          cParamSelTriCov         */
/*                                  */
/* ******************************** */

cParamSelTriCov::cParamSelTriCov() :
    mGreedy    (true),
    mNbMax     (8),
    mKappa     (0.3),
    mGainThr   (5e-3),
    mNbCells   (6),
    mBsHLim    (0.15),
    mAttenDens (3.0)
{
}

/* ******************************** */
/*                                  */
/*          cGridCoverTri           */
/*                                  */
/* ******************************** */

cGridCoverTri::cGridCoverTri(const std::vector<cPt2dr> & aVPts,int aNbCellsMax) :
    mP0   (aVPts.empty() ? cPt2dr(0,0) : aVPts.at(0)),
    mStep (0.0),
    mNbX  (1),
    mNbY  (1)
{
    MMVII_INTERNAL_ASSERT_always(!aVPts.empty(),"cGridCoverTri : no point");
    MMVII_INTERNAL_ASSERT_always(aNbCellsMax>0,"cGridCoverTri : NbCellsMax must be >0");

    cPt2dr aP1 = mP0;
    for (const auto & aPt : aVPts)
    {
        mP0 = cPt2dr(std::min(mP0.x(),aPt.x()),std::min(mP0.y(),aPt.y()));
        aP1 = cPt2dr(std::max(aP1.x(),aPt.x()),std::max(aP1.y(),aPt.y()));
    }
    cPt2dr aSz = aP1 - mP0;
    mStep = std::max(aSz.x(),aSz.y()) / aNbCellsMax;
    MMVII_INTERNAL_ASSERT_always(mStep>0,"cGridCoverTri : all points identical");

    // as MMV1 : longest side has NbCellsMax cells, the other proportionally less (at least 1)
    mNbX = std::min(aNbCellsMax,std::max(1,round_up(aSz.x()/mStep)));
    mNbY = std::min(aNbCellsMax,std::max(1,round_up(aSz.y()/mStep)));
}

size_t cGridCoverTri::NbCells() const
{
    return size_t(mNbX*mNbY);
}

size_t cGridCoverTri::IndexOf(const cPt2dr & aPt) const
{
    int anIx = std::min(mNbX-1,std::max(0,round_down((aPt.x()-mP0.x())/mStep)));
    int anIy = std::min(mNbY-1,std::max(0,round_down((aPt.y()-mP0.y())/mStep)));
    return size_t(anIy*mNbX + anIx);
}

std::vector<tREAL8> cGridCoverTri::Count(const std::vector<cPt2dr> & aVPts) const
{
    std::vector<tREAL8> aRes(NbCells(),0.0);
    for (const auto & aPt : aVPts)
        aRes.at(IndexOf(aPt)) += 1.0;
    return aRes;
}

/* ******************************** */
/*                                  */
/*      cSelecTripletsCoverage      */
/*                                  */
/* ******************************** */

cSelecTripletsCoverage::cSelecTripletsCoverage
    (
        const cParamSelTriCov &     aParam,
        const std::vector<cPt2dr> & aVPts2,
        const cPt3dr &              aCenterB,
        tREAL8                      aDepthH
        ) :
    mParam   (aParam),
    mGrid    (aVPts2,aParam.mNbCells),
    mNb2     (mGrid.Count(aVPts2)),
    mDensP   (mNb2.size(),0.0),
    mCenterB (aCenterB),
    mDepthH  (aDepthH)
{
    MMVII_INTERNAL_ASSERT_always(mDepthH>0,"cSelecTripletsCoverage : depth must be >0");

    tREAL8 aNbMax = *std::max_element(mNb2.begin(),mNb2.end()); // >0, aVPts2 not empty
    for (size_t aK=0 ; aK<mNb2.size() ; aK++)
        mDensP.at(aK) = std::sqrt(mNb2.at(aK)/aNbMax);
}

tREAL8 cSelecTripletsCoverage::CoverOfDensity(tREAL8 aDens,tREAL8 aAttenDens)
{
    return ((aAttenDens+1)/aAttenDens) * (aAttenDens*aDens) / (aAttenDens*aDens+1);
}

tREAL8 cSelecTripletsCoverage::GainOfBsH(tREAL8 aBsH,tREAL8 aBsHLim)
{
    return aBsH / (aBsH+aBsHLim);
}

tREAL8 cSelecTripletsCoverage::GainBsH(const cPt3dr & aC1,const cPt3dr & aC2) const
{
    return GainOfBsH(Norm2(aC1-aC2)/mDepthH,mParam.mBsHLim);
}

void cSelecTripletsCoverage::AddCandidate(const cCandTriCov & aCand)
{
    cCandInternal aCI;
    aCI.mId       = aCand.mId;
    aCI.mCenter   = aCand.mCenter;
    aCI.mSelected = false;

    // B/H gain of C : the worst of its two bases, with A (centre 0) and with B
    tREAL8 aGamma = std::min(GainBsH(cPt3dr(0,0,0),aCI.mCenter),GainBsH(mCenterB,aCI.mCenter));

    std::vector<tREAL8> aNb3 = mGrid.Count(aCand.mVPts3);
    for (size_t aK=0 ; aK<aNb3.size() ; aK++)
    {
        // triple points are not a subset of the pair points (MMV1 NB23), hence the min
        tREAL8 aDens = std::min(1.0,aNb3.at(aK)/std::max(1.0,mNb2.at(aK)));
        aCI.mCover.push_back(CoverOfDensity(aDens,mParam.mAttenDens));
        aCI.mGain.push_back(aCI.mCover.back()*aGamma);
    }
    mVCand.push_back(aCI);
}

tREAL8 cSelecTripletsCoverage::TotalGain(const cCandInternal & aCand) const
{
    tREAL8 aSum = 0.0;
    for (size_t aK=0 ; aK<aCand.mGain.size() ; aK++)
        aSum += aCand.mGain.at(aK) * mDensP.at(aK);
    return aSum / mGrid.NbCells();
}

void cSelecTripletsCoverage::UpdateGain(cCandInternal & aCand,const cCandInternal & aSel)
{
    tREAL8 aGammaCS = GainBsH(aCand.mCenter,aSel.mCenter);
    for (size_t aK=0 ; aK<aCand.mGain.size() ; aK++)
    {
        tREAL8 aDS  = aSel.mCover.at(aK);
        tREAL8 aCap = aCand.mCover.at(aK) * (aGammaCS*aDS + mParam.mKappa*(1-aDS));
        aCand.mGain.at(aK) = std::min(aCand.mGain.at(aK),aCap);
    }
}

std::vector<std::pair<int,tREAL8>> cSelecTripletsCoverage::Select()
{
    std::vector<std::pair<int,tREAL8>> aRes;

    // pass-through : every candidate, in insertion order, with its initial gain
    if (! mParam.mGreedy)
    {
        for (auto & aCand : mVCand)
        {
            aCand.mSelected = true;
            aRes.push_back(std::pair<int,tREAL8>(aCand.mId,TotalGain(aCand)));
        }
        return aRes;
    }

    while (aRes.size() < mParam.mNbMax)
    {
        // best remaining candidate ; strict > : first one in insertion order wins a tie
        cCandInternal * aBest = nullptr;
        tREAL8 aBestGain = -1.0;
        for (auto & aCand : mVCand)
        {
            if (! aCand.mSelected)
            {
                tREAL8 aGain = TotalGain(aCand);
                if (aGain > aBestGain)
                {
                    aBestGain = aGain;
                    aBest = & aCand;
                }
            }
        }
        if ((aBest==nullptr) || (aBestGain < mParam.mGainThr))
            break;

        aBest->mSelected = true;
        aRes.push_back(std::pair<int,tREAL8>(aBest->mId,aBestGain));
        for (auto & aCand : mVCand)
            if (! aCand.mSelected)
                UpdateGain(aCand,*aBest);
    }
    return aRes;
}

/* ******************************** */
/*                                  */
/*              Bench               */
/*                                  */
/* ******************************** */

namespace
{
typedef std::vector<std::pair<int,tREAL8>> tResSTC;

/// integer points (x,y), aX0<=x<aX1, 0<=y<aNbY ; with aOnlyEven, only those with x+y even
std::vector<cPt2dr> BenchSTC_Pts(int aX0,int aX1,int aNbY,bool aOnlyEven=false)
{
    std::vector<cPt2dr> aRes;
    for (int aX=aX0 ; aX<aX1 ; aX++)
        for (int aY=0 ; aY<aNbY ; aY++)
            if ((! aOnlyEven) || ((aX+aY)%2==0))
                aRes.push_back(cPt2dr(aX,aY));
    return aRes;
}

/// check ids in order and, for greedy runs, non-increasing gains
void BenchSTC_Check(const tResSTC & aRes,const std::vector<int> & aVIdExp,bool isGreedy,const std::string & aMes)
{
    MMVII_INTERNAL_ASSERT_bench(aRes.size()==aVIdExp.size(),aMes+" : bad number of selected");
    for (size_t aK=0 ; aK<aRes.size() ; aK++)
    {
        MMVII_INTERNAL_ASSERT_bench(aRes.at(aK).first==aVIdExp.at(aK),aMes+" : bad id/order");
        if (isGreedy && (aK>0))
        {
            MMVII_INTERNAL_ASSERT_bench(aRes.at(aK).second<=aRes.at(aK-1).second,aMes+" : gain increased");
        }
    }
}

/// run a selection on the 12x12 scene, A=0, B=(1,0,0), H=10
tResSTC BenchSTC_Run(const cParamSelTriCov & aParam,const std::vector<cCandTriCov> & aVCand)
{
    cSelecTripletsCoverage aSel(aParam,BenchSTC_Pts(0,12,12),cPt3dr(1,0,0),10.0);
    for (const auto & aCand : aVCand)
        aSel.AddCandidate(aCand);
    return aSel.Select();
}
}

void Bench_SelecTripletsCoverage(cParamExeBench & aParam)
{
    if (! aParam.NewBench("SelecTripletsCoverage")) return;

    const tREAL8 aEps = 1e-10;
    const std::vector<cPt2dr> aVAll   = BenchSTC_Pts(0,12,12);
    const std::vector<cPt2dr> aVLeft  = BenchSTC_Pts(0,6,12);
    const std::vector<cPt2dr> aVRight = BenchSTC_Pts(6,12,12);

    // ---- T1 : formulas and grid
    MMVII_INTERNAL_ASSERT_bench(std::abs(cSelecTripletsCoverage::CoverOfDensity(0.0,3.0))<aEps,"T1 : D(0)");
    MMVII_INTERNAL_ASSERT_bench(std::abs(cSelecTripletsCoverage::CoverOfDensity(1.0,3.0)-1.0)<aEps,"T1 : D(1)");
    MMVII_INTERNAL_ASSERT_bench(std::abs(cSelecTripletsCoverage::CoverOfDensity(0.5,3.0)-0.8)<aEps,"T1 : D(0.5)");
    MMVII_INTERNAL_ASSERT_bench(std::abs(cSelecTripletsCoverage::GainOfBsH(0.15,0.15)-0.5)<aEps,"T1 : gamma(Lim)");
    MMVII_INTERNAL_ASSERT_bench(std::abs(cSelecTripletsCoverage::GainOfBsH(0.0,0.15))<aEps,"T1 : gamma(0)");

    {
        cGridCoverTri aGrid(aVAll,6);
        MMVII_INTERNAL_ASSERT_bench(aGrid.NbCells()==36,"T1 : square grid size");
        for (const auto & aNb : aGrid.Count(aVAll))
        {
            MMVII_INTERNAL_ASSERT_bench(aNb==4,"T1 : square grid count");
        }
        MMVII_INTERNAL_ASSERT_bench(aGrid.IndexOf(cPt2dr(-5,-5))==0,"T1 : clamp low");
        MMVII_INTERNAL_ASSERT_bench(aGrid.IndexOf(cPt2dr(100,100))==35,"T1 : clamp high");

        cGridCoverTri aGridElong(BenchSTC_Pts(0,12,4),6);  // box 11x3 -> 6x2 cells
        MMVII_INTERNAL_ASSERT_bench(aGridElong.NbCells()==12,"T1 : elongated grid size");
    }

    cParamSelTriCov aPDef;
    cPt3dr aC3(0.5,1.5,0);
    tREAL8 aGamma3 = cSelecTripletsCoverage::GainOfBsH(Norm2(aC3)/10.0,0.15); // |C3|=|C3-B|
    {
        // full cover, P=1 everywhere : gain = gamma ; half density : gain = 0.8 gamma
        tResSTC aRes = BenchSTC_Run(aPDef,{{1,aC3,aVAll}});
        BenchSTC_Check(aRes,{1},true,"T1 full");
        MMVII_INTERNAL_ASSERT_bench(std::abs(aRes.at(0).second-aGamma3)<aEps,"T1 : full gain");

        aRes = BenchSTC_Run(aPDef,{{1,aC3,BenchSTC_Pts(0,12,12,true)}});
        BenchSTC_Check(aRes,{1},true,"T1 half");
        MMVII_INTERNAL_ASSERT_bench(std::abs(aRes.at(0).second-0.8*aGamma3)<aEps,"T1 : half gain");
    }

    // ---- T3a : C1 better than C2 initially, but redundant with C3 (close centre, same half)
    std::vector<cCandTriCov> aV3a {{1,cPt3dr(0.5,1.3,0),aVLeft},{2,cPt3dr(0.5,-1,0),aVRight},{3,aC3,aVLeft}};
    // ---- T3b : C1 far from C3 ; with Kappa=0.3, C1 beats C2 whose own B/H is capped
    std::vector<cCandTriCov> aV3b {{1,cPt3dr(0.5,0.5,0),aVLeft},{2,cPt3dr(0.5,-1,0),aVRight},{3,aC3,aVLeft}};

    cParamSelTriCov aPK1;
    aPK1.mKappa = 1.0;
    BenchSTC_Check(BenchSTC_Run(aPK1,aV3a),{3,2,1},true,"T3a K=1");
    BenchSTC_Check(BenchSTC_Run(aPDef,aV3a),{3,2,1},true,"T3a K=0.3");
    BenchSTC_Check(BenchSTC_Run(aPK1,aV3b),{3,2,1},true,"T3b K=1");
    BenchSTC_Check(BenchSTC_Run(aPDef,aV3b),{3,1,2},true,"T3b K=0.3");

    // ---- T4 : stop rules ; C4 has no triple point (gain 0)
    std::vector<cCandTriCov> aV4 = aV3b;
    aV4.push_back({4,cPt3dr(0.5,2,0),{}});
    {
        cParamSelTriCov aP = aPK1;
        aP.mNbMax = 2;
        BenchSTC_Check(BenchSTC_Run(aP,aV4),{3,2},true,"T4 NbMax");

        aP = aPK1;
        aP.mGainThr = 0.2;   // C1 gain 0.16 at third step
        BenchSTC_Check(BenchSTC_Run(aP,aV4),{3,2},true,"T4 GainThr");

        BenchSTC_Check(BenchSTC_Run(aPK1,aV4),{3,2,1},true,"T4 zero cover");
    }

    // ---- T5 : pass-through returns the whole pool, insertion order, initial gains
    {
        cParamSelTriCov aP;
        aP.mGreedy = false;
        tResSTC aRes = BenchSTC_Run(aP,aV4);
        BenchSTC_Check(aRes,{1,2,3,4},false,"T5");
        MMVII_INTERNAL_ASSERT_bench(std::abs(aRes.at(2).second-0.5*aGamma3)<aEps,"T5 : initial gain");
        MMVII_INTERNAL_ASSERT_bench(std::abs(aRes.at(3).second)<aEps,"T5 : zero gain");
    }

    aParam.EndBench();
}

} // MMVII