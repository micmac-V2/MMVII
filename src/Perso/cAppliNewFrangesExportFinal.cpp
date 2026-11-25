#include "cAppliNewFranges.h"



namespace MMVII
{

namespace NS_FrangesDetect
{

void cAppliNewFrange::OneIterCurvAndLabel(const tDIm& aDIm, tREAL8 aZoom)
{
    // initial attribution use  the labeling

    cPt2di aSzExport = aDIm.Sz();

    for (auto & aML : mVMaxLoc)
    {
        aML.mCurve = std::vector<tREAL8>(aSzExport.y(),-1);
        aML.mYMin =  1e8;
        aML.mYMax = -1e8;
        aML.mVLimX0 = std::vector<tREAL8>(aSzExport.y(),mSzRed.x());
        aML.mVLimX1 = std::vector<tREAL8>(aSzExport.y(),0);

    }

    for (int anY =0 ; anY<aSzExport.y() ; anY++)
    {
        std::vector<cWeightAv<tREAL8>> aVAvgX(mVMaxLoc.size());
        for (int anX=0 ; anX<aSzExport.x(); anX++)
        {
            cPt2di aPixFull(anX,anY);
            cPt2di aPixLab = mDImLabelEnd->Proj(ToI(ToR(aPixFull)/aZoom));

            int aLabel = mDImLabelEnd->GetV(aPixLab);
            if ( (aLabel>=0) && (aLabel<(int) mVMaxLoc.size()) )
            {
                aVAvgX.at(aLabel).Add(mDImBin->GetV(aPixFull),anX);
                if ((aLabel>=0) && (aLabel<(int)mVMaxLoc.size()))
                {
                   UpdateMin(mVMaxLoc.at(aLabel).mYMin,aPixFull.y());
                   UpdateMax(mVMaxLoc.at(aLabel).mYMax,aPixFull.y());

                   UpdateMin(mVMaxLoc.at(aLabel).mVLimX0.at(anY),(tREAL8) aPixFull.x());
                   UpdateMax(mVMaxLoc.at(aLabel).mVLimX1.at(anY),(tREAL8) aPixFull.x());
                }
            }
        }

        for (size_t aKC=0 ; aKC<mVMaxLoc.size() ; aKC++)
        {
            if (aVAvgX.at(aKC).SW())
            {
               mVMaxLoc.at(aKC).mCurve.at(anY) = aVAvgX.at(aKC).Average();
            }
        }
    }
}

void cAppliNewFrange::OneIterRelabel()
{

    if (mVMaxLoc.size()<=1)
        return;

    mDImLabelEnd->InitCste(mVMaxLoc.size());


    for (size_t aLabel=0 ; aLabel<mVMaxLoc.size() ; aLabel++)
    {
        cMaxLocSum & aML = mVMaxLoc.at(aLabel);
        if (aML.mYMin<=aML.mYMax)
        {
            for (int anY =aML.mYMin ; anY<=aML.mYMax ; anY++)
            {
                tREAL8 aXC =  mVMaxLoc.at(aLabel).mCurve.at(anY);
                tREAL8 aX0 ;//= aXC;
                tREAL8 aX1 ;//= aXC;

                if (aLabel==0)
                {
                   const cMaxLocSum & aNext = mVMaxLoc.at(aLabel+1);
                   if (aNext.YIsInit(anY))
                   {
                      aX0 = aXC + (aXC-aNext.mCurve.at(anY)) / 2.0;
                   }
                   else
                   {
                       aX0 = aXC / 2.0;
                   }
                }
                else
                {
                     cMaxLocSum & aPrec = mVMaxLoc.at(aLabel-1);
                     if (aPrec.YIsInit(anY))
                         aX0 = (aPrec.mCurve.at(anY) +aXC) / 2.0;
                     else
                         aX0 = aXC/2.0;
                }


                if (aLabel==mVMaxLoc.size()-1)
                {
                    const cMaxLocSum & aPrec = mVMaxLoc.at(aLabel-1);
                    if (aPrec.YIsInit(anY))
                    {
                       aX1 = aXC + (aXC-aPrec.mCurve.at(anY)) / 2.0;
                    }
                    else
                    {
                        aX1 = (aXC+mSzRed.x()) /2.0;
                    }
                }
                else
                {
                    cMaxLocSum & aNext = mVMaxLoc.at(aLabel+1);
                    if (aNext.YIsInit(anY))
                        aX1 = (aNext.mCurve.at(anY) +aXC) / 2.0;
                    else
                        aX1 = (mSzRed.x() +aXC) / 2.0;
                }

                int iX0 = std::clamp(round_ni(aX0),0,mSzRed.x()-1);
                int iX1 = std::clamp(round_ni(aX1),0,mSzRed.x()-1);

                for (int anX=iX0; anX<=iX1 ; anX++)
                {
                     mDImLabelEnd->SetV(cPt2di(anX,anY),aLabel);
                }

                 aML.mVLimX0.at(anY) = iX0;
                 aML.mVLimX1.at(anY) = iX1;

            }
        }
    }
}


void cAppliNewFrange::FixPbLastFrange()
{
    //  it possible that a small part of franges on the right has not created a line
    for (int anY=0 ; anY<mSzRed.y() ; anY++)
    {
        int anX0 = mVMaxLoc.back().mX1+ round_ni(mDImIntegr->GetV(anY)*2.0);

        for (int anX= anX0 ; anX<mSzRed.x() ; anX++)
        {
            mDImLabelEnd->SetV(cPt2di(anX,anY),mVMaxLoc.size());
        }
    }
}


void cAppliNewFrange::MakeFinalAttribution(const tDIm& aDIm, tREAL8 aZoom)
{

    mImLabelEnd  = cIm2D<tU_INT1>(mSzRed);
    mDImLabelEnd =  &mImLabelEnd.DIm();
    mDImLabelEnd->InitCste(255);


    // initial attribution use  the labeling
    for (const auto aPix :  *mDImRegion)
    {
        int aReg = mDImRegion->GetV(aPix);
        if (aReg>=0)
        {
            tVertexReg & aVert1 = mGrReg.VertexOfNum(aReg);
            int aNumCC = aVert1.Attr().mNumCC0;
            if (aNumCC<0)
                aNumCC = 128;
            mDImLabelEnd->SetV(aPix,aNumCC);
            if ((aNumCC>=0) && (aNumCC<(int)mVMaxLoc.size()))
            {
               UpdateMin(mVMaxLoc.at(aNumCC).mYMin,aPix.y());
               UpdateMax(mVMaxLoc.at(aNumCC).mYMax,aPix.y());
            }
        }
    }
    cPt2di aSzExport = aDIm.Sz();

    //  it possible that a small part of franges on the right has not created a line


    FixPbLastFrange();


    mImBin = cIm2D<tREAL4>(aSzExport,nullptr,eModeInitImage::eMIA_Null);
    mDImBin = & mImBin.DIm();
    mDImBin->InitCste(128);

    for (int anY =0 ; anY<aSzExport.y() ; anY++)
    {
        cComputeStdDev<tREAL8> aStat;
        for (int anX=0 ; anX<aSzExport.x(); anX++)
        {
            cPt2di aPix(anX,anY);
            aStat.Add(aDIm.GetV(aPix));
        }

        aStat = aStat.Normalize(1e-5);

       // int aYTgt = anY/ aZoom;
      //  aYTgt = std::clamp(aYTgt,0,mDImTgt->Sz()-1);
      //  tREAL8 aTgt = std::abs(mDImTgt->GetV(aYTgt));
        //tREAL8 aMulErf = 1.0 +  4.0 / (1.0+aTgt);
        tREAL8 aMulErf = 4.0;
        tREAL8 aLowT = 0.0;// -1/3.0;
        tREAL8 aPowLow = 1.0;

        for (int anX=0 ; anX<aSzExport.x(); anX++)
        {
            cPt2di aPix(anX,anY);
            tREAL8 aV0 = aDIm.GetV(aPix);
            tREAL8 aV = aStat.NormalizedVal(aV0);

            aV = std::erf(aV*aMulErf);
            if (aV<aLowT)
            {
                aV = aLowT - std::pow(aLowT-aV,aPowLow);
            }
            //aV = aStat.GaussianNormalizedRank(aV,4.0);

            mDImBin->SetVTrunc(aPix,127*(1+aV));
        }
    }

   // mDImBin->ToFile("TotoBin.tif");

    OneIterCurvAndLabel(aDIm,aZoom);
   // OneIterCurvAndLabel();
   // OneIterCurvAndLabel();

}


};  // NS_FrangesDetect
}; //  namespace MMVII
