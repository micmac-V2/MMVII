//#include "MMVII_PCSens.h"

// Test commit

#include "cAppliNewFranges.h"




namespace MMVII
{

namespace NS_FrangesDetect
{


/* ******************************************** */
/*                                              */
/*             cConnComp                        */
/*                                              */
/* ******************************************** */

cConnComp::cConnComp(std::vector<cPt2di> & aVPts,const tSeg2dr & aSeg,cPt2dr aTg1,cPt2dr aTg2,bool isHor,tREAL8 aScore,bool isOk) :
    mPts   (aVPts),
    mSeg   (aSeg),
    mTg1   (aTg1),
    mTg2   (aTg2),
    mIsHor (isHor),
    mScore (aScore),
    mIsOk  (isOk)
{
}


/* ******************************************** */
/*                                              */
/*             cAppliNewFrange                  */
/*                                              */
/* ******************************************** */

/*
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
*/

cPt2di cAppliNewFrange::NewBackTrackExtre(cPt2di  aP0,const std::vector<cPt2di> & aVPts,tREAL8 aRay) const
{
    tREAL8 aR2 = Square(aRay);
    cWhichMin<cPt2di,tREAL8> aWMin;

    for (const auto & aPt : aVPts)
        aWMin.Add(aPt,std::abs(SqN2(aP0-aPt)-aR2));

    return aWMin.IndexExtre();
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


     cPt2dr aP1 =  ToR(NewBackTrackExtre(ToI(aSeg.P1()),aVCC,10.0));
     cPt2dr aP2 =  ToR(NewBackTrackExtre(ToI(aSeg.P2()),aVCC,10.0));

   //  StdOut() << "SEG=" << aP1-aSeg.P1() << aP2-aSeg.P2() << "\n";

     if ( (Norm2(aP1-aSeg.P1())<3.0) || (Norm2(aP2-aSeg.P2())<3.0))
         return false;


     aP1 = VUnit(aSeg.P1() - aP1);
     aP2 = VUnit(aSeg.P2() - aP2);

    // StdOut() << "P1P2=" << aP1 << aP2 << "\n";


     cConnComp * aCC = new cConnComp(aVCC,aSeg,aP1,aP2,isHoriz,aScore,isOk);
   //  cConnComp(const tSeg2dr & aSeg,bool isHor) :

     mVectCC.push_back(aCC);

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


void cAppliNewFrange::ComplemCCHoriz()
{
    tREAL8 aPropLow = 0.6;
    tREAL8 aThreshold = ((aPropLow*mRadiomBackGround) + (1-aPropLow) * mDImRadFrange->GetV(mYC));
    int aNbMaxLoc = 40;
    tREAL8 aDistCC = 100.0;
    //int aNbX = mSzRed.x();


    for (int anX=aNbMaxLoc ; anX< mSzRed.x()-1-aNbMaxLoc ; anX++)
    {
        cPt2di aPt(anX,mYC);
        tREAL8 aVC = mDImRedBlur->GetV(aPt) ;
        bool isMaxLoc = aVC > aThreshold;

        if (isMaxLoc)
        {
            for (int aDx=1 ; (aDx<aNbMaxLoc) && isMaxLoc ; aDx++)
            {
                tREAL8 aVPDx =  mDImRedBlur->GetV(cPt2di(anX+aDx,mYC));
                tREAL8 aVPMx =  mDImRedBlur->GetV(cPt2di(anX-aDx,mYC));

                 isMaxLoc = isMaxLoc && (aVC>=aVPMx) && (aVC>aVPDx);
            }
        }

        if (isMaxLoc)
        {
           for (size_t aKC=0 ; (aKC<mVectCC.size()) && isMaxLoc ; aKC++)
           {
               const cConnComp & aCC = *mVectCC.at(aKC);
               if (aCC.mIsHor && aCC.mIsOk)
               {
                   cSegment2DCompiled<tREAL8> aSeg(aCC.mSeg);
                   tREAL8 aD = aSeg.DistClosedSeg(ToR(aPt));
                   isMaxLoc = isMaxLoc && (aD>aDistCC);
               }
           }
        }


      //(std::vector<cPt2di> & aVPts,const tSeg2dr & aSeg,cPt2dr aTg1,cPt2dr aTg2,bool isHor,tREAL8 aScore,bool isOk) :

        if (isMaxLoc)
        {
            int aBorder = std::min(mSzRed.y()/4,aNbMaxLoc);
            int anY0 = mYC;
            while ((anY0>aBorder) &&  (mDImRedBlur->GetV(cPt2di(anX,anY0))>aThreshold) )
                anY0--;

            int anY1 = mYC;
            while (anY1<mSzRed.y() - aBorder &&  (mDImRedBlur->GetV(cPt2di(anX,anY1))>aThreshold) )
                anY1++;

            std::vector<cPt2di>  aVPt;
            for (int anY=anY0 ; anY<=anY1 ; anY++)
                aVPt.push_back(cPt2di(anX,anY));

            tSeg2dr aSeg(cPt2dr(anX,anY0),cPt2dr(anX,anY1));
            cConnComp* aCC = new cConnComp(aVPt,aSeg,cPt2dr(0,-1),cPt2dr(0,1),true,1.0,true);
            mVectCC.push_back(aCC);
        }

    }

   // for (int a)
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

    ComplemCCHoriz();
}



};  // namespace NS_FrangesDetect
}; //  namespace MMVII




