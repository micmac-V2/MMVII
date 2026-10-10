//#include "MMVII_PCSens.h"

// Test commit

#include "cAppliNewFranges.h"




namespace MMVII
{

namespace NS_FrangesDetect
{

cPt2di cAppliNewFrange::OneStepBackTrackExtre(const cPt2di & aP0,const cPt2di & aPLast,const cDataIm2D<tU_INT1> & aImMax) const
{
    cPt2di aRes = aPLast;
    const std::vector<cPt2di> & aVec8Neigh =  Alloc8Neighbourhood();

    int aNb=0;
    int aCoul = aImMax.GetV(aPLast);

    for (const auto & aNeigh : aVec8Neigh)
    {
        cPt2di aPt = aP0 + aNeigh;
        if ((aImMax.GetV(aPt)==aCoul) && (aPt!=aPLast))
        {
            aRes = aPt;
            aNb++;
        }
    }
    if (aNb!=1) aRes = aPLast;

    return aRes;
}



cPt2di cAppliNewFrange::BackTrackExtre(cPt2di  aP0,const cDataIm2D<tU_INT1> & aImMax) const
{
   cPt2di aP1 =    OneStepBackTrackExtre(aP0,aP0,aImMax);

   for (int aK=0 ; (aK< 10) && (aP1!=aP0) ; aK++ )
   {
       cPt2di aP2 =  OneStepBackTrackExtre(aP1,aP0,aImMax);
       aP0 = aP1;
       aP1 = aP2;
   }

   return aP1;
}


bool cAppliNewFrange::NewCC(bool isHoriz,const cPt2di& aPix)
{
    const std::vector<cPt2di> & a8Neigh =  Alloc8Neighbourhood();
    cDataIm2D<tU_INT1> & aImMax = isHoriz ? *mDImMaxHor : * mDImMaxVert;
    std::vector<cPt2di> aVCC;

    ConnectedComponent (aVCC,aImMax ,a8Neigh, aPix,1,2);

    bool isInMax = false;
    if (mHasMask && isHoriz)
    {
       for (const auto & aPix : aVCC)
       {
          if (mImMask.DIm().GetV(aPix))
          {
              isInMax = true;
          }
       }
       FakeUseIt(isInMax);
    }

    tREAL8 aThrs = isHoriz ? 5 : 20 ;
    if (aVCC.size()<aThrs)
        return false;

    //cBox2di aBox(cTplBox<int,2>::FromVect(aVCC));
    cTplBox<int,2>  aBox =cTplBox<int,2>::FromVect(aVCC,true);

    bool aY0Inf =  (aBox.P0().y() <= mYC);
    bool aY1Sup =  (aBox.P1().y() > mYC);
    // MMVII_INTERNAL_ASSERT_always(aY0Inf<=aY1Sup," Ordeeerr in box");

     bool doCross = (aY0Inf && aY1Sup);


// FakeUseIt(doCross);
     cPt2di aDir(0,1);
     if (isHoriz)
     {
         // if dont cross
        if (!doCross)
            return false;
     }
     else
     {
      //   if (doCross)
       //     return false;
         aDir = (aBox.P0().y()<mYC)  ? cPt2di(-1,0) : cPt2di(1,0);
     }


     cWhichMinMax<cPt2di,tREAL8>  aWMM;
     cWeightAv<tREAL8,tREAL8> aWScore;

     for (const auto & aPt : aVCC)
     {
         aWMM.Add(aPt,Scal(aDir,aPt));

         tREAL8 aV = mDImRedBlur->GetV(aPt);
         tREAL8 aAmpl = mDImRadFrange->GetV(aPt.y()) -mRadiomBackGround;
         aV = (aV-mRadiomBackGround)/aAmpl;
         aV = std::min(1.0,std::max(aV,0.0));

         aWScore.Add(1.0,aV);
     }

     tREAL8 aScore = aWScore.Average();
     tREAL8 aScoreMixte = aScore * aVCC.size();
     bool isOk =    (aScore > 0.5)
                 || (aScoreMixte > 15)
                 || ((!isHoriz) && (aVCC.size() > 150))
            ;


     tSeg2dr aSeg(ToR(aWMM.IndMin()),ToR(aWMM.IndMax()));


     cPt2dr aP1 =  ToR(BackTrackExtre(ToI(aSeg.P1()),aImMax));
     cPt2dr aP2 =  ToR(BackTrackExtre(ToI(aSeg.P2()),aImMax));

   //  StdOut() << "SEG=" << aP1-aSeg.P1() << aP2-aSeg.P2() << "\n";

     if ( (Norm2(aP1-aSeg.P1())<3.0) || (Norm2(aP2-aSeg.P2())<3.0))
         return false;


     aP1 = VUnit(aSeg.P1() - aP1);
     aP2 = VUnit(aSeg.P2() - aP2);

    // StdOut() << "P1P2=" << aP1 << aP2 << "\n";


     cConnComp aCC(aVCC,aSeg,aP1,aP2,isHoriz,aScore,isOk);
   //  cConnComp(const tSeg2dr & aSeg,bool isHor) :

     mListCC.push_back(aCC);

     return true;
}

void cAppliNewFrange::ConnecCompMaxLoc(bool isHoriz)
{
    cDataIm2D<tU_INT1> & aImMax = isHoriz ? *mDImMaxHor : * mDImMaxVert;

    for (const auto aPix : aImMax)
    {
        if (aImMax.GetV(aPix)==1)
        {
            NewCC(isHoriz,aPix);
        }
    }
}


bool cAppliNewFrange::IsMax(cPt2di aP0,cPt2di aDp0,int aNb) const
{
    auto * anDIm = mDImRedBlur;

    tREAL4 aV0 = anDIm->GetV(aP0);
    for (int aK=1 ; aK<aNb ; aK++)
    {
        cPt2di aDP = aDp0 * aK;

        if (     ( anDIm->GetV(aP0 + aDP) > aV0)
              || ( anDIm->GetV(aP0 - aDP) > aV0)
           )
        {
            return false;
        }
    }

    return true;
}

void cAppliNewFrange::MakeImageMaxLoc()
{
    mImMaxHor = cIm2D<tU_INT1>(mSzRed,nullptr,eModeInitImage::eMIA_Null);
    mDImMaxHor = & (mImMaxHor.DIm());
    mImMaxVert = cIm2D<tU_INT1>(mSzRed,nullptr,eModeInitImage::eMIA_Null);
    mDImMaxVert = & (mImMaxVert.DIm());


    int aNbMaxX=8;
    int aNbMaxY=3;

    for (const auto aPix : mDImRed->Interior(1+std::max(aNbMaxX,aNbMaxY)))
    {
         if ( IsMax(aPix,cPt2di(1,0),aNbMaxX))
        {
             mDImMaxHor->SetV(aPix,1);
        }

        if (  IsMax(aPix,cPt2di(0,1),aNbMaxY))
        {
              mDImMaxVert->SetV(aPix,1);
        }
    }

    ConnecCompMaxLoc(true);
    ConnecCompMaxLoc(false);

}



};  // namespace NS_FrangesDetect
}; //  namespace MMVII




