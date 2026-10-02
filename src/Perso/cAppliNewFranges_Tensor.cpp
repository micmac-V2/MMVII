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
tREAL8 cAppliNewFrange::QualityCenterAntiSymetry(int aY0,int aSzWindow)
{
    tREAL8 aSumSign = 0;
    tREAL8 aSumAbs = 0;

    for (int aDY=1 ; aDY<=aSzWindow; aDY++)
    {
        tREAL8 aV1 = mDImTeta->GetV(aY0-aDY);
        tREAL8 aV2 = mDImTeta->GetV(aY0+aDY);

        tREAL8 aW = aSzWindow - std::abs(aDY);

        aSumSign += Square(aV1+aV2) * aW;
        aSumAbs += (Square(aV1) + Square(aV2)) *aW;
    }

    return aSumSign / aSumAbs;
}

/* Compute the center as the value minimzing QualityCenterAntiSymetry */

int cAppliNewFrange::ComputeCenterByAntiSym()
{
    // Compute size of window, adapt it height of image
    int aSzY = round_up(100 / mZoomRed);
    aSzY = std::min(aSzY,mSzRed.y()/4);

    cWhichMin<int,tREAL8> aMinSym;
    int aEndY =  mDImTeta->Sz()-(aSzY+1);

    for (int aY0=aSzY ; aY0 <aEndY ; aY0++)
        aMinSym.Add(aY0,QualityCenterAntiSymetry(aY0,aSzY));

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



void cAppliNewFrange::DoTensorProcessing()
{
    // Compute the deriche gradient
    cImGrad<tElIm>  aGrad = Deriche(*mDImRed,mDerFactZ1*mZoomRed);
    tDIm & aDGx = *(aGrad.mDGx);
    tDIm & aDGy = *(aGrad.mDGy);

    // Resize the images of tensor-2D and
    mImTx =  tIm(mSzRed);
    mDImTx = &(mImTx.DIm());
    mImTy  = tIm(mSzRed);
    mDImTy = &(mImTy.DIm());
    mImMaxLoc = cIm2D<tU_INT1> (mSzRed,nullptr,eModeInitImage::eMIA_Null);
    mDImMaxLoc = &(mImMaxLoc.DIm());

    // Images of tensor accumulation
     cIm1D<tREAL8> aPop(mSzRed.y(),nullptr,eModeInitImage::eMIA_Null);
     cIm1D<tREAL8> aSumTx(mSzRed.y(),nullptr,eModeInitImage::eMIA_Null);
     cIm1D<tREAL8> aSumTy(mSzRed.y(),nullptr,eModeInitImage::eMIA_Null);


     // Parse all point to store 2D-tensor and accumlate in 1D images
     for (const auto & aPix : aDGx)
     {
         cPt2dr aGrad(aDGx.GetV(aPix),aDGy.GetV(aPix));
         cPt2dr aRhoTeta = ToPolar(aGrad,0.0);

         tREAL8 aRho = aRhoTeta.x();
         tREAL8 aTeta = aRhoTeta.y();
         cPt2dr aTens = FromPolar(aRho,2.0*aTeta);

         mDImTx->SetV(aPix,aTens.x());
         mDImTy->SetV(aPix,aTens.y());

         tREAL8 aWeight = 1.0;

          aPop.DIm().AddV(aPix.y(),aRho*aWeight);
          aSumTx.DIm().AddV(aPix.y(),aTens.x()*aWeight);
          aSumTy.DIm().AddV(aPix.y(),aTens.y()*aWeight);
     }

     // Regularize 2D images of tensor
     ExpFilterOfStdDev(*mDImTx,5,3.0);
     ExpFilterOfStdDev(*mDImTy,5,3.0);

     // tentative extraction of centers of lines as maxima of grey
     // in direction of tensor
     int aBorder = 5;
     std::vector<cPt2dr> aVN;
     for (const auto & aPix : mDImMaxLoc->Interior(aBorder))
     {
         cPt2dr aTens (mDImTx->GetV(aPix),mDImTy->GetV(aPix));
         cPt2dr aRhoTeta = ToPolar(aTens,0.0);
         cPt2dr aDirTens = FromPolar(1.0,aRhoTeta.y()/2.0);
         tREAL8 aV0 =mDImRedBlur->GetV(aPix);
         if (      (aV0>mDImRedBlur->GetVBL(ToR(aPix)+aDirTens))
               &&   (aV0>mDImRedBlur->GetVBL(ToR(aPix)-aDirTens))
            )
         {
             mDImMaxLoc->SetV(aPix,1);
         }
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

     for(const auto anY : aPop.DIm())
     {
         cPt2dr aTens(aSumTx.DIm().GetV(anY),aSumTy.DIm().GetV(anY));
         cPt2dr aRhoTeta = ToPolar(aTens,0.0);
         tREAL8 aTeta = aRhoTeta.y() /2.0;
         tREAL8 aTgt = tan(aTeta);

         mDImTgt->SetV(anY,aTgt);
         mDImTeta->SetV(anY,aTeta);

     }

     // --------  Extract center ------------------------
     mYC = ComputeCenterByAntiSym();


     // --------------- Make integrale images ---------------------
     mImIntegr  = cIm1D<tREAL8>(mSzRed.y(),nullptr,eModeInitImage::eMIA_Null);
     mDImIntegr = &(mImIntegr.DIm());

     OneWayIntegrateTangent(-1,0);
     OneWayIntegrateTangent(1,mSzRed.y());

}


};  // NS_FrangesDetect
}; //  namespace MMVII

