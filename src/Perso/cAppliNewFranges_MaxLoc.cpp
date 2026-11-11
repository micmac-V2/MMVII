//#include "MMVII_PCSens.h"

// Test commit

#include "cAppliNewFranges.h"


namespace MMVII
{

namespace NS_FrangesDetect
{


/* ******************************************** */
/*                                              */
/*               cCCMax                         */
/*                                              */
/* ******************************************** */

cAVertReg::cAVertReg(const  std::vector<cPt2di> & aVPts,int aNum) :
    mVPts (aVPts),
    mNum  (aNum),
    mCoul ((53*mNum)%256,(137*(mNum+1))%256,(109*(mNum+2))%256)
{
}

cAOrientReg::cAOrientReg()
{
}

cASymReg::cASymReg() :
    mMaxVal (-1e10)
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

    size_t aThBig = (aLabel==eTypePix::eMinTens) ? 100 : 40;
    size_t aThSmall= (aLabel==eTypePix::eMinTens) ? 50 : 10;
    tREAL8 aThAvg = (aLabel==eTypePix::eMinTens) ? 0.7  : 0.5;
    tREAL8 aExpAvg = (aLabel==eTypePix::eMinTens) ? 0.25  : 1.0;

    if (aVCC.size()> aThBig)
        return aLabel;

    if (aVCC.size()< aThSmall)
        return eTypePix::eTypeSmalCC;

    cWeightAv<tREAL8,tREAL8> aWScore;
    for (const auto & aPt : aVCC)
    {
        tREAL8 aV = mDImRedBlur->GetV(aPt);
        tREAL8 aAmpl = mDImRadFrange->GetV(aPt.y()) -mRadiomBackGround;
        aV = (aV-mRadiomBackGround)/aAmpl;
        aV = std::min(1.0,std::max(aV,0.0));
        aV = std::pow(aV,aExpAvg);
        aWScore.Add(1.0,aV);
    }
    tREAL8 aScore = std::pow(aWScore.Average(),1/aExpAvg);
    if (aLabel==eTypePix::eMinTens)
        aScore = 1.0-aScore;

    if (aScore< aThAvg)
        return  eTypePix::eTypeBadAvgCC;


    return aLabel;

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


#endif
}

void cAppliNewFrange::ConnecCompMaxLoc( eTypePix aLabelIn)
{
    // Create the vectext corresponding to exterioir
    if (aLabelIn== eTypePix::eMaxTens)
    {
        std::vector<cPt2di> aVCC;
        mVirtExtVert =   mGrReg.NewVertex(cAVertReg(aVCC,mGrReg.NbVertex()));
    }


    for (const auto aPix : *mDImLabel)
    {
        if (LabelIsOn(aPix,aLabelIn) && LabelTmpIsOff(aPix))
        {
            std::vector<cPt2di> aVCC;
            eTypePix aLabelOut =   NewCC(aLabelIn,aVCC,aPix);

            if (aLabelIn==aLabelOut) // The CC was "successful"
            {

                if (aLabelIn==eTypePix::eMaxTens)
                {
                    mGrReg.NewVertex(cAVertReg(aVCC,mGrReg.NbVertex()));
                }
            }
            else // fail, reset value of point to label representing the reason of fail
            {
                for (const auto & aPt : aVCC)
                    SetLabel(aPt,aLabelOut);
            }
        }
    }

    // The Tmp Flag must be reset for next connect comp
    for (const auto aPix : *mDImLabel)
        SetLabelTmpOff(aPix);
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


void cAppliNewFrange::AnalyseCC()
{

     ConnecCompMaxLoc(eTypePix::eMaxTens);
     ConnecCompMaxLoc(eTypePix::eMinTens);


}



};  // namespace NS_FrangesDetect
}; //  namespace MMVII




