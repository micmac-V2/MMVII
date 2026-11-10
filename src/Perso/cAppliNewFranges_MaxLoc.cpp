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


/**
 * For extracting direction extremity, conpute the point with is at a given distance
 * (or rathet the closest at this distance)
 */

cPt2di cAppliNewFrange::NewBackTrackExtre(cPt2di  aP0,const std::vector<cPt2di> & aVPts,tREAL8 aRay) const
{
    tREAL8 aR2 = Square(aRay);
    cWhichMin<cPt2di,tREAL8> aWMin;

    for (const auto & aPt : aVPts)
        aWMin.Add(aPt,std::abs(SqN2(aP0-aPt)-aR2));

    return aWMin.IndexExtre();
}


eTypePix cAppliNewFrange::NewCC(eTypePix aLabel,std::vector<cPt2di>  & aVCC,const cPt2di& aPix)
{
    const std::vector<cPt2di> & a8Neigh =  Alloc8Neighbourhood();
    std::vector<cPt2di>  aSeed{aPix};

    FlagBitConnectedComponent (aVCC,*mDImLabel ,a8Neigh, aSeed,Lab2I(aLabel),LabTmp2I());

    /// Test is in mask
 /*   bool isInMax = false;
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
    }*/



#if (0)

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
  //   tREAL8 aScoreMixte = aScore * aVCC.size();
    /* bool isOk =    (aScore > 0.5)
                 || (aScoreMixte > 15)
                 || ((!isHoriz) && (aVCC.size() > 150))*/

     bool isOk =  (aVCC.size() > 100) ;
     if (!isOk)
         return false;


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

#endif
     return aLabel;
}

void cAppliNewFrange::ConnecCompMaxLoc( eTypePix aLabelIn)
{
    //cDataIm2D<tU_INT1> & aImMax = isHoriz ? *mDImMaxHor : * mDImMaxVert;

    for (const auto aPix : *mDImLabel)
    {
        if (LabelIsOn(aPix,aLabelIn) && LabelTmpIsOff(aPix))
        {
            std::vector<cPt2di> aVCC;
            eTypePix aLabelOut =   NewCC(aLabelIn,aVCC,aPix);
            if (aLabelIn==aLabelOut)
            {
            }
            else
            {
                for (const auto & aPt : aVCC)
                    SetLabelOn(aPt,aLabelOut);
            }
        }
    }

}


bool cAppliNewFrange::AddFor8Cx(const cDataIm2D<tU_INT1> & aImMax,const cPt2di& aPix) const
{
    static const std::vector<cPt2di> & aV4 = Alloc4Neighbourhood();

    if (aImMax.GetV(aPix)==2)
        return false;

    int aNbN = 0;
    for (const auto & aNeigh : aV4)
      if (aImMax.GetV(aPix+aNeigh)==2)
           aNbN++;

    return aNbN==2;
}
void cAppliNewFrange::Make8Cx(cDataIm2D<tU_INT1> & aImMax)
{
    std::vector<cPt2di> aVPt;
    for (const auto & aPix : aImMax.Interior(1))
    {
        if (AddFor8Cx(aImMax,aPix))
            aVPt.push_back(aPix);
    }

    //StdOut() <<  "Make8CxMake8CxMake8Cx " << aVPt.size() << "\n";
    for (const auto & aPix : aVPt)
        aImMax.SetV(aPix,2);
}


bool cAppliNewFrange::IsMinMax(cPt2di aP0,cPt2di aDp0,int aNb,bool isMax) const
{
    int aSign = isMax ? 1 : -1;
    auto * anDIm = mDImRedBlur;

    tREAL4 aV0 = anDIm->GetV(aP0);
    for (int aK=1 ; aK<aNb ; aK++)
    {
        cPt2di aDP = aDp0 * aK;

        if (     ( aSign*anDIm->GetV(aP0 + aDP) > aSign*aV0)
              || ( aSign*anDIm->GetV(aP0 - aDP) >= aSign*aV0)
           )
        {
            return false;
        }
    }

    return true;
}

void cAppliNewFrange::MakeImageMaxLoc()
{
 /*   mImMaxHor = cIm2D<tU_INT1>(mSzRed,nullptr,eModeInitImage::eMIA_Null);
    mDImMaxHor = & (mImMaxHor.DIm());
    mImMaxVert = cIm2D<tU_INT1>(mSzRed,nullptr,eModeInitImage::eMIA_Null);
    mDImMaxVert = & (mImMaxVert.DIm());*/


    int mDistVert2C = 25;
    int aNbMaxX=8;
    int aNbMaxY=3;

    for (const auto aPix : mDImRed->Interior(1+std::max(aNbMaxX,aNbMaxY)))
    {
         if ( IsMinMax(aPix,cPt2di(0,1),aNbMaxY,false)  && (std::abs(aPix.y()-mYC)>=mDistVert2C))
        {         // cIm2D<tU_INT1> mImMaxLocTD;           ///<  local maxima in tensor direction
             // cDataIm2D<tU_INT1>* mDImMaxLocTD;     ///< data

            // mDImMaxHor->SetV(aPix,1);
             SetLabelOn(aPix,eTypePix::eMinVert);
        }

        if (  IsMinMax(aPix,cPt2di(0,1),aNbMaxY,true) && (std::abs(aPix.y()-mYC)>=mDistVert2C))
        {
           //   mDImMaxVert->SetV(aPix,1);
              SetLabelOn(aPix,eTypePix::eMaxVert);
        }
    }


     ConnecCompMaxLoc(eTypePix::eMinVert);
     ConnecCompMaxLoc(eTypePix::eMaxVert);
     ConnecCompMaxLoc(eTypePix::eMaxTens);




  //  Make8Cx(*mDImMaxHor);

    //ConnecCompMaxLoc(*mDImMaxVert,false);

    //ComplemCCHoriz();
}



};  // namespace NS_FrangesDetect
}; //  namespace MMVII




