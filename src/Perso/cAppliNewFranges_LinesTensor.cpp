//#include "MMVII_PCSens.h"

// Test commit

#include "cAppliNewFranges.h"




namespace MMVII
{


namespace NS_FrangesDetect
{
#if (0)
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



tREAL8 AdaptTetaPos(tREAL8 aTeta,tREAL8 aMaxTeta)
{
    if (aTeta < - M_PI/4)
        aTeta += M_PI;

    return std::min(aTeta,aMaxTeta);
//    return std::max(0.0,std::min(aTeta,aMaxTeta));
}

tREAL8 AdaptTetaNeg(tREAL8 aTeta,tREAL8 aMaxTeta)
{
    return - AdaptTetaPos(-aTeta,aMaxTeta);
}
#endif
bool cAppliNewFrange::IsSaddlePoint(const cPt2di & aPix) const
{
    if (mDImRedBlur->Insideness(aPix)<=1)
        return false;
   return  FlagSup8Neigh(*mDImRedBlur,aPix).NbConComp() >=4;

}

cPt2dr cAppliNewFrange::TensorImage(cPt2di aPix) const
{
   return ToR(cTensor<tElIm>::Tens2Grad(mImTens.Grad(aPix),1.0));

 //  return  FromPolar(1.0,mDImTeta->GetV(aPix.y()));


}


cPt2dr cAppliNewFrange::NextPt(const cPt2dr & aP1,const cPt2dr& aP2) const
{
    static std::vector<cPt2dr> aVDir;
    if (aVDir.empty())
    {
        tREAL8 aStep = 0.15;
        tREAL8 aTetaMax = 0.3;
        int aNbTeta = round_up(aTetaMax/aStep);
        aStep = aTetaMax / aNbTeta;
        for (int aKTeta=-aNbTeta; aKTeta<=aNbTeta ; aKTeta++)
        {
            aVDir.push_back(FromPolar(1.0,aStep*aKTeta));
        }
    }
    cWhichMax<cPt2dr,tREAL8> aMaxP;

    cPt2dr aDirNext = VUnit(aP2-aP1);
    for (const auto&  aDir : aVDir)
    {
        cPt2dr aNewP = aP2+aDir * aDirNext;
        tREAL8 aVal = mDImRedBlur->DefGetVBL(aNewP,-1e8);
        aMaxP.Add(aNewP,aVal);
    }

    return aMaxP.IndexExtre();
}


void cAppliNewFrange:: ProlongateExtre(std::vector<cPt2dr>  & aRes, cPt2di aPix,eTypePix aLabel) const
{
    aRes.clear();
    if (! LabelIsOn(aPix,aLabel))
        return;

    cFlag8Neigh aFlag8N;
    for (size_t aKN=0 ; aKN<8 ; aKN++)
    {
          if (LabelIsOn(aPix+FreemanV8[aKN],aLabel))
              aFlag8N.AddNeigh(aKN);
    }

    if (aFlag8N.NbConComp()!=2)
        return;

    tREAL8 aVal = mDImRedBlur->GetV(aPix);
    tREAL8 aThrs = aVal + (aVal-mRadiomBackGround);


    aRes.push_back(ToR(aPix));
    cPt2dr aDirT = TensorImage(aPix);

    cWhichMax<cPt2dr,tREAL8> aClosTens;
    for (const auto& aDN8 : Alloc8Neighbourhood())
    {
        cPt2di aN8 = aPix + aDN8;
        if ( LabelIsOn(aN8,aLabel))
        {
            aClosTens.Add(ToR(aDN8),std::abs(Scal(aDirT,VUnit(ToR(aDN8)))));
        }
    }
    cPt2dr aDelta =  aClosTens.IndexExtre();
    cPt2dr aNext = ToR(aPix) - aDelta;
    aRes.push_back(aNext);

    for (int aK=0 ; aK<6 ; aK++)
    {
        cPt2dr aNew = NextPt(aRes.at(aRes.size()-2),aRes.at(aRes.size()-1));
        aRes.push_back(aNew);

        cPt2di aINew = ToI(aNew);
        if ((mDImLabel->Insideness(aINew)<3) ||  (mDImRedBlur->GetV(aINew) > aThrs))
        {
            aRes.clear();
            return ;
        }



        for (const auto & aNeigh : Alloc8Neighbourhood())
        {
             if (LabelIsOn(aINew+aNeigh,aLabel))
                 return ;
        }
    }

    aRes.clear();

}



bool cAppliNewFrange::IsMaxLocTens(const cPt2di & aPix,const cPt2dr & aDir) const
{
    tREAL8 aV0 = mDImRedBlur->GetV(aPix);
    tREAL8 aVMin1 = 1e10;
    tREAL8 aVMin2 = 1e10;


    int aKMinGlob = 15;
    int aKMinLoc =   5;


    for (int aK=1 ; aK<aKMinGlob ; aK++)
    {
        cPt2dr aP1 = ToR(aPix) + aDir * (tREAL8)aK;
        if (mDImRedBlur->InsideBL(aP1))
        {
            tREAL8 aV1 =  mDImRedBlur->GetVBL(aP1);
            if  ((aV1>aV0) && (aK<aKMinLoc))
                return false;
            UpdateMin(aVMin1,aV1);
        }

        cPt2dr aP2 = ToR(aPix) - aDir * (tREAL8)aK;
        if (mDImRedBlur->InsideBL(aP2))
        {
            tREAL8 aV2 =  mDImRedBlur->GetVBL(aP2);
            if ((aV2>=aV0) &&  (aK<aKMinLoc))
                return false;
            UpdateMin(aVMin2,aV2);
        }
    }

    tREAL8 aContrastMin = (mDImRadFrange->GetV(aPix.y()) - mRadiomBackGround) * 0.05;

    if ( ((aV0-aVMin1)<aContrastMin) || ( (aV0-aVMin2)<aContrastMin) )
        return false;

    return true;
}

#if (0)
void cAppliNewFrange::DoTensorProcessing()
{
    // Compute the deriche gradient
   // mImGrad = Deriche(*mDImRed,mDerFactZ1*mZoomRed);

    tIm aImBlurGrad = mImRed.Dup();
    ExpFilterOfStdDev(aImBlurGrad.DIm(),5,mSigTensors.mSigmaGradX,mSigTensors.mSigmaGradY);
    mImGrad = Sobel(aImBlurGrad.DIm());

    mImTens = cTensor<tElIm>::Grad2Tens(mImGrad);
    // Regularize 2D images of tensor
    ExpFilterOfStdDev(*mImTens.mDGx,5,mSigTensors.mSigmaTensX,mSigTensors.mSigmaTensY);
    ExpFilterOfStdDev(*mImTens.mDGy,5,mSigTensors.mSigmaTensX,mSigTensors.mSigmaTensY);



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

     //  first compute on tensor to avoid the discontinuity arrond pi/2
     mYC = ComputeCenterByAntiSym(true,aVDirASym,aVWeight);

     //  Adapt teta to avoid discont (
     if (1)
     {
         tREAL8 aMaxTan = 20.0;
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
    /* tREAL8 aS2T = Square(mSigTensors.mSigmaTeta);
      StdOut() << " FACT SIG TETA " <<  FactExpFromSigma2(aS2T)
               << " " << FactExpFromSigma2(aS2T) / aS2T  << "\n";*/


     if (0)
     {
      /*   tIm1D aImSigma(mSzRed.y());

         for ( int anY=0 ; anY<mSzRed.y() ; anY++)
         {
             int aNSig = 10;
             int aY0 =
         }*/
     }
     else
         ExpFilterOfStdDev(*mDImTeta,3,mSigTensors.mSigmaTeta);

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

     if (1)
     {
        ComputeMaxRadCenter();
     }
     else
     {
         ComputeSumOnIntegralCurve();
     }
     ExtractMaxSumIC();


     ComputeRadiomCste();
}
#endif

void cAppliNewFrange::ExtractTensorLines()
{

     // tentative extraction of centers of lines as maxima of grey
     // in direction of tensor
     int aBorder = 5;
     //std::vector<cPt2dr> aVN;
     std::vector<cPt2di> aVPixMax;

     for (const auto & aPix : mDImLabel->Interior(aBorder))
     {
         cPt2dr aDirTens = ToR(cTensor<tElIm>::Tens2Grad(mImTens.Grad(aPix),1.0));



         aDirTens = FromPolar(1.0,mDImTeta->GetV(aPix.y()));

         if (IsMaxLocTens(aPix,aDirTens))
         {
            aVPixMax.push_back(aPix);
            SetLabelOn(aPix,eTypePix::eMaxTens);
          }
     }

     for (const auto & aPixMax : aVPixMax)
     {
         std::vector<cPt2dr> aVProl;
         ProlongateExtre(aVProl,aPixMax,eTypePix::eMaxTens);
         for (const auto & aProl : aVProl)
         {
             SetLabelOn(ToI(aProl),eTypePix::eMaxTens);
             SetLabelOn(ToI(aProl),eTypePix::eProlong);
         }
     }
}




};  // NS_FrangesDetect
}; //  namespace MMVII

