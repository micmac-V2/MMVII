//#include "MMVII_PCSens.h"

// Test commit

#include "cAppliNewFranges.h"




namespace MMVII
{


namespace NS_FrangesDetect
{

bool cAppliNewFrange::IsMaxLocSumIC(int anX) const
{
    int mBorder = 10;
    int mSzMax = 50;

    if ((anX<mBorder) || (anX>=mSzRed.x()-mBorder))
        return false;

    tREAL8 aVC = mDImProjInt->GetV(anX);
    tREAL8 aVMin = 1e10;

    for (int aDX=1 ; aDX<mSzMax ; aDX++)
    {
        int anX0 = anX -aDX;
        if (anX0>=0)
        {
            tREAL8 aV0 = mDImProjInt->GetV(anX0) ;
            if (aV0 >= aVC)
               return false;
            UpdateMin(aVMin,aV0);
        }

        int anX1 = anX + aDX;    
        if (anX1<mSzRed.x()-1)
        {
             tREAL8 aV1 = mDImProjInt->GetV(anX1) ;
             if (aV1>aVC)
                return false;
             UpdateMin(aVMin,aV1);
        }
    }

    // StdOut() << (aVC-aVMin) / mStdDevProjI << "\n";

    if ((aVC-aVMin) < 0.2 * mStdDevProjI)
        return false;
    return true;
}

void cAppliNewFrange::ExtractMaxSumIC()
{
    for (int anX=0 ; anX<mSzRed.x() ; anX++)
        if (IsMaxLocSumIC(anX))
        {
            cMaxLocSum aMLS;
            aMLS.mXC = anX;
            mVMaxLoc.push_back(aMLS);
        }

    std::vector<tREAL8> aVInterv;
    for (size_t aK1=1 ; aK1<mVMaxLoc.size(); aK1++)
    {
        tREAL8 aXLim = (mVMaxLoc.at(aK1-1).mXC + mVMaxLoc.at(aK1).mXC ) / 2.0;
        mVMaxLoc.at(aK1-1).mX1 = aXLim;
        mVMaxLoc.at(aK1).mX0 = aXLim;
        aVInterv.push_back(mVMaxLoc.at(aK1).mXC - mVMaxLoc.at(aK1-1).mXC);
    }

    mVMaxLoc.front().mX0 = 2*mVMaxLoc.front().mXC - mVMaxLoc.front().mX1;
    mVMaxLoc.back().mX1  = 2*mVMaxLoc.back().mXC  - mVMaxLoc.back().mX0 ;

    mIntervMax = NonConstMediane(aVInterv);
}



void cAppliNewFrange::ComputeSumOnIntegralCurve()
{

    /*  To take into account the delocalsation of the curve :
     *     - make a "strong" smoothing along line of ImRed save it in mDImSmoothTgt,
     *     - the smmothing depend of the line : higer when Y is far from YC
     *     - mehod :
     *         o parse lines
     *         o make a copy of each line, make the gaussian filter and save back
     */

    mImSmoothTgt = tIm(mSzRed);
    mDImSmoothTgt = &mImSmoothTgt.DIm();

    cIm1D<tREAL4> anImL(mSzRed.x());
    for (int anY=0 ; anY<mSzRed.y() ; anY++)
    {
        mDImRed->WriteLine(anImL.DIm(),anY);
        /*for (int anX=0 ; anX<mSzRed.x() ; anX++)
        {
            anImL.DIm().SetV(anX,mDImRed->GetV(cPt2di(anX,anY)));
        }*/

        tREAL8 aSig0 = 3.0;
        tREAL8 aSigMin = 30.0;
        tREAL8 aSigm = std::sqrt(Square(aSig0)+ Square(mDImTgt->GetV(anY)));
        aSigm *= aSigMin/aSig0;
        ExpFilterOfStdDev(anImL.DIm(),2,aSigm);

        mDImSmoothTgt->ReadLine(anImL.DIm(),anY);
/*
        for (int anX=0 ; anX<mSzRed.x() ; anX++)
        {
            mDImSmoothTgt->SetV(cPt2di(anX,anY),anImL.DIm().GetV(anX));
        }*/
    }


    tREAL8 aDistW = 50.0;
    mImProjInt = tIm1D(mSzRed.x());
    mDImProjInt = &mImProjInt.DIm();

    cIm1D<tU_INT2> aImCpt(mSzRed.x(),nullptr,eModeInitImage::eMIA_Null);
    int aThrsNb = mSzRed.y()/10;
    mMinDPI = 1e10;
    mMaxDPI = 0;
    for (int anX=0 ; anX<mSzRed.x() ; anX++)
    {
        cWeightAv<tREAL8> aAvg;
        int aNb =0;
        for (int anY=0 ; anY<mSzRed.y() ; anY++)
        {
            tREAL8 aXTgt = anX+ mDImIntegr->GetV(anY);
            cPt2dr aPt(aXTgt,anY);
            if (mDImSmoothTgt->InsideBL(aPt))
            {
                aNb++;
                tREAL8 aW =  aDistW / (aDistW + std::abs(anY-mYC));
                aAvg.Add(aW,mDImSmoothTgt->GetVBL(aPt));
            }
        }
        if (aNb>aThrsNb)
        {
            mDImProjInt->SetV(anX,aAvg.Average());
            UpdateMin(mMinDPI,aAvg.Average());
            UpdateMax(mMaxDPI,aAvg.Average());

        }
        aImCpt.DIm().SetV(anX,aNb);
    }

    for (int anX=0 ; anX<mSzRed.x() ; anX++)
    {
        if (aImCpt.DIm().GetV(anX)<=aThrsNb)
        {
            mDImProjInt->SetV(anX,mMinDPI);
        }
    }
}

void      cAppliNewFrange::ComputeMaxRadCenter()
{
    mImProjInt = tIm1D(mSzRed.x());
    mDImProjInt = &mImProjInt.DIm();

    tIm aImFilter = mImRed.Dup();
    ExpFilterOfStdDev(aImFilter.DIm(),3,30.0,10.0);

    cStdStatRes aStat;
    mMinDPI = 1e10;
    mMaxDPI = 0;

    for (const auto & aP : *mDImProjInt)
    {
        tREAL8 aVal = aImFilter.DIm().GetV(cPt2di(aP.x(),mYC));
        mDImProjInt->SetV(aP,aVal);
        aStat.Add(aVal);
        UpdateMin(mMinDPI,aVal);
        UpdateMax(mMaxDPI,aVal);
    }
    mMinDPI = aStat.Min();
    mMaxDPI = aStat.Max();
    mStdDevProjI = aStat.DevStd();
}



};  // NS_FrangesDetect
}; //  namespace MMVII

