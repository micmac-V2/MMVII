//#include "MMVII_PCSens.h"

// Test commit

#include "cAppliNewFranges.h"




namespace MMVII
{


namespace NS_FrangesDetect
{

/*   Compute how good Y0 is a center, quantifying the antisymetry of teta,
 *   the formula is :
 *                    Sum (  (Teta(y)+Teta(-y))^2 )  / Sum( Teta(y)^2 +Teta(-y)^2)
 */
tREAL8 cAppliNewFrange::QualityCenterAntiSymetry(bool First,int aY0,int aSzWindow,const std::vector<cPt2dr>&aVDir,const std::vector<tREAL8>& aVW)
{

    cWeightAv<tREAL8,tREAL8> aAvgSign;
    cWeightAv<tREAL8,tREAL8> aAvg2;

    for (int aDY=1 ; aDY<=aSzWindow; aDY++)
    {
        int aY1 = aY0-aDY;
        int aY2 = aY0+aDY;
        tREAL8 aWeight =  (aSzWindow - std::abs(aDY));


        if (First)
        {
           cPt2dr aP1 = aVDir.at(aY1);
           cPt2dr aP2 = aVDir.at(aY2);
           aAvgSign.Add(aWeight,SqN2(aP1+aP2));
           aAvg2.Add(aWeight,SqN2(aP1)+SqN2(aP2));
        }
        else
        {
            tREAL8 aTeta1 = mDImTeta->GetV(aY1);
            tREAL8 aTeta2 = mDImTeta->GetV(aY2);

             aAvgSign.Add(aWeight,Square(aTeta1+aTeta2));
             aAvg2.Add(aWeight,Square(aTeta1)+Square(aTeta2));

        }

    }

  /*  if (First)
        return aAvgSign.Average();
    else */
        return  aAvgSign.Average() / aAvg2.Average();
}

/* Compute the center as the value minimzing QualityCenterAntiSymetry */

int cAppliNewFrange::ComputeCenterByAntiSym(bool First,const std::vector<cPt2dr>& aVDir,const std::vector<tREAL8>& aVW)
{



    // Compute size of window, adapt it height of image
    int aSzY = round_up(100 / mZoomRed);
    aSzY = std::min(aSzY,mSzRed.y()/4);

    cWhichMin<int,tREAL8> aMinSym;
    int aEndY =  mDImTeta->Sz()-(aSzY+1);

    for (int aY0=aSzY ; aY0 <aEndY ; aY0++)
        aMinSym.Add(aY0,QualityCenterAntiSymetry(First,aY0,aSzY,aVDir,aVW));

    return aMinSym.IndexExtre();
}


void cAppliNewFrange::OneWayIntegrateTangent(int aDy,int aYLim)
{
    for (int aY=mYC+aDy  ; aY!= aYLim ; aY+= aDy)
    {
        tREAL8 aPreVal =  mDImIntegr->GetV(aY-aDy);
        tREAL8 aPreTgt =  mDImTgt->GetV(aY-aDy);
        tREAL8 aCurTgt =  mDImTgt->GetV(aY);

        mDImIntegr->SetV(aY,aPreVal+(-aDy)*(aPreTgt+aCurTgt)/2.0);
    }
}

void cAppliNewFrange::ComputeProjIntegrX()
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
    for (int anY=0 ; anY<mSzRed.y() ; anY++)
    {
        tIm1D anImL(mSzRed.x());
        for (int anX=0 ; anX<mSzRed.x() ; anX++)
        {
            anImL.DIm().SetV(anX,mDImRed->GetV(cPt2di(anX,anY)));
        }

        tREAL8 aSig0 = 3.0;
        tREAL8 aSigMin = 30.0;
        tREAL8 aSigm = std::sqrt(Square(aSig0)+ Square(mDImTgt->GetV(anY)));
        aSigm *= aSigMin/aSig0;
        ExpFilterOfStdDev(anImL.DIm(),2,aSigm);

        for (int anX=0 ; anX<mSzRed.x() ; anX++)
        {
            mDImSmoothTgt->SetV(cPt2di(anX,anY),anImL.DIm().GetV(anX));
        }
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

tREAL8 AdaptTetaPos(tREAL8 aTeta,tREAL8 aMaxTeta)
{
    if (aTeta < - M_PI/4)
        aTeta += M_PI;

    return std::max(0.0,std::min(aTeta,aMaxTeta));
}

tREAL8 AdaptTetaNeg(tREAL8 aTeta,tREAL8 aMaxTeta)
{
    return - AdaptTetaPos(-aTeta,aMaxTeta);
}





void cAppliNewFrange::DoTensorProcessing()
{
    // Compute the deriche gradient
   // mImGrad = Deriche(*mDImRed,mDerFactZ1*mZoomRed);

    tIm aImBlurGrad = mImRed.Dup();
    ExpFilterOfStdDev(aImBlurGrad.DIm(),5,2.0,2.0);
    mImGrad = Sobel(aImBlurGrad.DIm());

    mImTens = cTensor<tElIm>::Grad2Tens(mImGrad);
    // Regularize 2D images of tensor
    ExpFilterOfStdDev(*mImTens.mDGx,5,20,5);
    ExpFilterOfStdDev(*mImTens.mDGy,5,20,5);



    mImCumNormT = tIm1D(mSzRed.y(),nullptr,eModeInitImage::eMIA_Null);
    mDImCumNormT = & mImCumNormT.DIm();

    // Images of tensor accumulation
     cIm1D<tREAL8> aPop(mSzRed.y(),nullptr,eModeInitImage::eMIA_Null);
     cIm1D<tREAL8> aSumTx(mSzRed.y(),nullptr,eModeInitImage::eMIA_Null);
     cIm1D<tREAL8> aSumTy(mSzRed.y(),nullptr,eModeInitImage::eMIA_Null);


     // Parse all point to store 2D-tensor and accumlate in 1D images
     for (const auto & aPix : *mDImRed)
     {
         cPt2dr   aTens = mImTens.GradR(aPix);//Tensor<tREAL8>::Grad2Tens(mImGrad.GradR(aPix));
         tREAL8 aWeight = 1.0;

         aPop.DIm().AddV(aPix.y(),1.0);
         aSumTx.DIm().AddV(aPix.y(),aTens.x()*aWeight);
         aSumTy.DIm().AddV(aPix.y(),aTens.y()*aWeight);
         mDImCumNormT->AddV(aPix.y(),SqN2(aTens));
     }



    // ------------ Filter 1D accumlator, dans make average -------------

     ExpFilterOfStdDev(aPop.DIm()  ,5,mSigmaTensZ1/mZoomRed);
     ExpFilterOfStdDev(aSumTx.DIm(),5,mSigmaTensZ1/mZoomRed);
     ExpFilterOfStdDev(aSumTy.DIm(),5,mSigmaTensZ1/mZoomRed);

     DivImageInPlace(aSumTx.DIm(),aSumTx.DIm(),aPop.DIm());  // divide
     DivImageInPlace(aSumTy.DIm(),aSumTy.DIm(),aPop.DIm());

     // ------- Transformate the 1D-tensor  accum in teta/tan

     mImTgt = cIm1D<tREAL8>(mSzRed.y(),nullptr,eModeInitImage::eMIA_Null);
     mDImTgt = & (mImTgt.DIm());
     mImTeta =  cIm1D<tREAL8>(mSzRed.y(),nullptr,eModeInitImage::eMIA_Null);
     mDImTeta = &(mImTeta.DIm());

     std::vector<cPt2dr> aVDirASym;
     std::vector<tREAL8>  aVWeight;

     for(const auto anY : aPop.DIm())
     {
         cPt2dr aTens(aSumTx.DIm().GetV(anY),aSumTy.DIm().GetV(anY));
         cPt2dr aRhoTeta = ToPolar(aTens,0.0);
         tREAL8 aRho = aRhoTeta.x();
         tREAL8 aTeta = aRhoTeta.y();
         aVDirASym.push_back(FromPolar(1.0,aTeta) - cPt2dr(1.0,0)) ;
         aVWeight.push_back(aRho);

         mDImTeta->SetV(anY,aTeta/2.0);
     }


     // --------  Extract center ------------------------
     mYC = ComputeCenterByAntiSym(true,aVDirASym,aVWeight);

     if (1)
     {
         tREAL8 aMaxTan = 15.0;
         tREAL8 aMaxTeta = atan(aMaxTan);
         for (const auto anY : aPop.DIm())
         {
             tREAL8 aTeta = mDImTeta->GetV(anY);
             aTeta = (anY.x()<= mYC) ? AdaptTetaPos(aTeta,aMaxTeta) : AdaptTetaNeg(aTeta,aMaxTeta);
             mDImTeta->SetV(anY,aTeta);
         }
     }

     //tREAL8 aSigma = std::min(20.0,);
     mYC = ComputeCenterByAntiSym(false,aVDirASym,aVWeight);

    // tREAL8 aSigma = std::min(15.0,mSzRed.y()/50.0);
     ExpFilterOfStdDev(*mDImTeta,3,15.0);

     for(const auto anY : aPop.DIm())
     {
         mDImTgt->SetV(anY,tan( mDImTeta->GetV(anY)));
     }
   //  mYC = ComputeCenterByAntiSym();

     // --------------- Make integrale images ---------------------
     mImIntegr  = cIm1D<tREAL8>(mSzRed.y(),nullptr,eModeInitImage::eMIA_Null);
     mDImIntegr = &(mImIntegr.DIm());

     OneWayIntegrateTangent(-1,-1);
     OneWayIntegrateTangent(1,mSzRed.y());

     ComputeProjIntegrX();


     {
         // Resize ImMaxLoc
         mImMaxLocTD = cIm2D<tU_INT1> (mSzRed,nullptr,eModeInitImage::eMIA_Null);
         mDImMaxLocTD = &(mImMaxLocTD.DIm());


         // tentative extraction of centers of lines as maxima of grey
         // in direction of tensor
         int aBorder = 5;
         std::vector<cPt2dr> aVN;
         for (const auto & aPix : mDImMaxLocTD->Interior(aBorder))
         {
             cPt2dr aDirTens = ToR(cTensor<tElIm>::Tens2Grad(mImTens.Grad(aPix),1.0));

             tREAL8 aV0 =mDImRedBlur->GetV(aPix);
             if (      (aV0>mDImRedBlur->GetVBL(ToR(aPix)+aDirTens))
                   &&   (aV0>mDImRedBlur->GetVBL(ToR(aPix)-aDirTens))
                )
             {
                 mDImMaxLocTD->SetV(aPix,1);
             }
         }
     }

}




};  // NS_FrangesDetect
}; //  namespace MMVII

