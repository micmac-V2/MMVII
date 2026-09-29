//#include "MMVII_PCSens.h"

#include "MMVII_Image2D.h"
#include "MMVII_ImageMorphoMath.h"
#include "MMVII_Sensor.h"
#include "MMVII_Ptxd.h"
#include "MMVII_2Include_Serial_Tpl.h"
#include "MMVII_Tpl_ElemStrToVal.h"
#include "MMVII_Tpl_Images.h"
#include "MMVII_Tpl_ElemFilterLocImages.h"

#include "MMVII_Tpl_GraphStruct.h"
#include "MMVII_Tpl_Graph_SubGraph.h"
#include "MMVII_Tpl_GraphAlgo_SPCC.h"




namespace MMVII
{

namespace NS_FrangesDetect
{

/**  An application for  extaction curves on interference images
 *   rather very specific...
 */
typedef tREAL4            tElIm;
typedef cIm2D<tElIm >     tIm;
typedef cDataIm2D<tElIm > tDIm;
typedef cIm1D<tREAL8>     tIm1D;
typedef cDataIm1D<tREAL8>     tDIm1D;


/* =============================================== */
/*                                                 */
/*                 cAppliNewFrange                   */
/*                                                 */
/* =============================================== */


class cAppliNewFrange : public cMMVII_Appli
{
     public :
        cAppliNewFrange(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec);
     private :
        //================= typedef part ===============================


        typedef cIm2D<tU_INT1 >    tImLabel;
        typedef cDataIm2D<tU_INT1> tDImLabel;

        //================================================================
        //       METHODS DECLARATION
        //================================================================

            //------------------------- overidding cMMVII_Appli ----------------
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;
         std::vector<cOneHelpSampleCmp>  Samples() const override;
        virtual ~cAppliNewFrange();

     private :

        void  DoOneImage(const std::string & aNameIm);
        std::string NameVisu(const std::string & aPref) const;

        void MakeImSimul();

        /** Compute images of local maxims, this will reduce "drastically" the number of potential point
         in shortest path approach (and allow more complexe lins with "big" jumps ) */
        void MakeImageMaxLoc();

        void MakeImTgt();

        tREAL8 CostSymTgt(int aY,int aSzY);
        int DetectSymByTgt();
        void DoIntegrale(int aDy,int aYLim);

        cPhotogrammetricProject  mPhProj;

        // ----------- Mandatory Args -----------
        std::string              mPatImage; /// Pattern of all images
        std::string              mNameIm;

        tREAL8    mZoomRed;     ///< Zoom for initial reduction
        tREAL8    mDerFactZ1;   ///< Factor of deriche gradient for zoom 1
        tREAL8    mSigmaTensZ1; ///< Sigma filter on tensor cumulated for zoom 1
        bool      mDoSimul;
        bool      mDoVisu;

        tIm      mImZ1;
        tDIm*    mDImZ1;

        bool IsMax(cPt2di aP0,cPt2di aDp,int aNb) const;


        //-------- These values are related to reduced image ----------------
        tIm      mImRed;
        tDIm*    mDImRed;
        cPt2di   mSzRed;
        tIm      mImRedBlur;
        tDIm*    mDImRedBlur;



         cIm2D<tU_INT1> mImMax;
         cDataIm2D<tU_INT1>* mDImMax;



        tIm1D    mImTgt;
        tDIm1D*  mDImTgt;
        int      mYC;
        tIm1D    mImIntegr;
        tDIm1D*  mDImIntegr;
        cRGBImage mImVisu;

};

    /* =================================================== */
    /*      Overiding of  cMMVII_Appli                     */
    /* =================================================== */


cAppliNewFrange::cAppliNewFrange(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
    cMMVII_Appli      (aVArgs,aSpec),
    mPhProj           (*this),
    mZoomRed          (2.0),
    mDerFactZ1        (0.5),
    mSigmaTensZ1      (10.0),
    mDoSimul          (false),
    mDoVisu           (true),
    mImZ1             (cPt2di(1,1)),
    mDImZ1            (nullptr),
    mImRed            (cPt2di(1,1)),
    mDImRed           (nullptr),
    mImRedBlur        (cPt2di(1,1)),
    mDImRedBlur       (nullptr),
    mImMax            (cPt2di(1,1)),
    mDImMax           (nullptr),
    mImTgt            (1),
    mDImTgt           (nullptr),
    mYC               (-1),
    mImIntegr         (1),
    mDImIntegr        (nullptr),
    mImVisu           (cPt2di(1,1))
{
}

cAppliNewFrange::~cAppliNewFrange()
{
}

cCollecSpecArg2007 & cAppliNewFrange::ArgObl(cCollecSpecArg2007 & anArgObl)
{
      return    anArgObl
            <<  Arg2007(mPatImage,"Name of input Image", {eTA2007::FileDirProj,{eTA2007::MPatFile,"0"}})
      ;
}

cCollecSpecArg2007 & cAppliNewFrange::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    anArgOpt
              << AOpt2007(mZoomRed,"Zoom","Zoom Red process init" ,{eTA2007::HDV})
              << AOpt2007(mDerFactZ1,"DericheF","Factor of deriche filter (for Zoom=1)" ,{eTA2007::HDV})
              << AOpt2007(mSigmaTensZ1,"SigmaTens","Sigma for avaragin tensor (Z=1)" ,{eTA2007::HDV})
              << AOpt2007(mDoSimul,"DoSimul","Make simulation/syntheic images" ,{eTA2007::HDV})

                 //  << AOpt2007(mSigCurv,"SigCurv","Sima for smoothig curve",{eTA2007::HDV})
            //  << AOpt2007(mIntY,"IntY","Interval for Y",{eTA2007::HDV})
            //  << AOpt2007(mNbIter,"NbIt","Number of iter initial smoothing",{eTA2007::HDV})
            //  << AOpt2007(mDoVisu,"DoVisu","Generate Visualisation ?",{eTA2007::HDV})
            //  << AOpt2007(mWidhMinAll,"MinWidth","Minimal witdh, general case",{eTA2007::HDV})
            //  << AOpt2007(mWidhMinBorder,"BorderMinWidth","Minimal witdh for border",{eTA2007::HDV})
            //  << AOpt2007(mMinHeightBorderRight,"MinHeightBorderRight","Minimal Heitgh for right border",{eTA2007::HDV})
            //  << AOpt2007(mParamStack,"Stack","Stacking parm [Pat,Nb,Mode] ",{{eTA2007::ISizeV,"[3,3]"}})
    ;

     return anArgOpt;
}


std::vector<cOneHelpSampleCmp>  cAppliNewFrange::Samples() const
{
   return
   {
       {"MMVII ExtractFranges Retiga_000000105.tif Sigma=[10,2] DoVisu=1"},
       {"MMVII ExtractFranges Retiga_.*.tif"}
   };
}

int cAppliNewFrange::Exe()
{
    mPhProj.FinishInit();

    if (RunMultiSet(0,0)) // Case several images in //
    {
       return ResultMultiSet();
    }
    DoOneImage(UniqueStr(0)); // Case 1 image (may be recalled from multiple)
    return EXIT_SUCCESS;
}

/* =================================================== */
/*              Specific functions                     */
/* =================================================== */
void cAppliNewFrange::MakeImSimul()
{
    tREAL8 aMulY=10.0;
    tREAL8 aPer = 300.0 / mZoomRed;
    tREAL8 aExp = 4.0;
    tREAL8 aMilY = mSzRed.y() / 2.0;
    for (const auto & aPix : *mDImRed)
    {
        tREAL8 aCorrY = (aPix.y()-aMilY) / aMilY;

        tREAL8 aPhase = aPix.x() -Square(aCorrY) * aMilY * aMulY;


        aPhase /= (aPer /(2*M_PI));
        tREAL8 aAmpl = std::max(0.0,(1+sin(aPhase)) /2.0);
        aAmpl = std::pow(aAmpl,aExp);

        aAmpl = aAmpl / (1 + std::pow(std::abs(aCorrY),2)*3.0 );

        mDImRed->SetV(aPix,aAmpl*255.0);
    }
}


tREAL8 cAppliNewFrange::CostSymTgt(int aY0,int aSzY)
{
    tREAL8 aSumSign = 0;
    tREAL8 aSumAbs = 0;

    for (int aDY=1 ; aDY<=aSzY; aDY++)
    {
        tREAL8 aV1 = mDImTgt->GetV(aY0-aDY);
        tREAL8 aV2 = mDImTgt->GetV(aY0+aDY);

        tREAL8 aW = aSzY - std::abs(aDY);

        aSumSign += Square(aV1+aV2) * aW;
        aSumAbs += (Square(aV1) + Square(aV2)) *aW;
    }

    return aSumSign / aSumAbs;
}

int cAppliNewFrange::DetectSymByTgt()
{
    int aSzY = round_up(100 / mZoomRed);
    aSzY = std::min(aSzY,mSzRed.y()/4);

    cWhichMin<int,tREAL8> aMinSym;

    int aEndY =  mDImTgt->Sz()-(aSzY+1);

    for (int aY0=aSzY ; aY0 <aEndY ; aY0++)
        aMinSym.Add(aY0,CostSymTgt(aY0,aSzY));

    return aMinSym.IndexExtre();
}

void cAppliNewFrange::DoIntegrale(int aDy,int aYLim)
{
    for (int aY=mYC+aDy  ; aY!= aYLim ; aY+= aDy)
    {
        tREAL8 aPreVal =  mDImIntegr->GetV(aY-aDy);
        tREAL8 aPreTgt =  mDImTgt->GetV(aY-aDy);
        tREAL8 aCurTgt =  mDImTgt->GetV(aY);

        mDImIntegr->SetV(aY,aPreVal+(-aDy)*(aPreTgt+aCurTgt)/2.0);
    }
}

bool cAppliNewFrange::IsMax(cPt2di aP0,cPt2di aDp0,int aNb) const
{
    tREAL4 aV0 = mDImRedBlur->GetV(aP0);
    for (int aK=1 ; aK<aNb ; aK++)
    {
        cPt2di aDP = aDp0 * aK;

        if (     ( mDImRedBlur->GetV(aP0 + aDP) > aV0)
              || ( mDImRedBlur->GetV(aP0 - aDP) > aV0)
           )
        {
            return false;
        }
    }

    return true;
}


void cAppliNewFrange::MakeImTgt()
{
    cImGrad<tElIm>  aGrad = Deriche(*mDImRed,mDerFactZ1*mZoomRed);
    tDIm & aDGx = *(aGrad.mDGx);
    tDIm & aDGy = *(aGrad.mDGy);

    // aDGx.ToFile("GradX.tif");

     cIm1D<tREAL8> aPop(mSzRed.y(),nullptr,eModeInitImage::eMIA_Null);
     cIm1D<tREAL8> aSumTx(mSzRed.y(),nullptr,eModeInitImage::eMIA_Null);
     cIm1D<tREAL8> aSumTy(mSzRed.y(),nullptr,eModeInitImage::eMIA_Null);


     for (const auto & aPix : aDGx)
     {
         cPt2dr aGrad(aDGx.GetV(aPix),aDGy.GetV(aPix));
         cPt2dr aRhoTeta = ToPolar(aGrad,0.0);

         tREAL8 aRho = aRhoTeta.x();
         tREAL8 aTeta = aRhoTeta.y();
         cPt2dr aTens = FromPolar(aRho,2.0*aTeta);

         tREAL8 aWeight = 1.0;

          aPop.DIm().AddV(aPix.y(),aRho*aWeight);
          aSumTx.DIm().AddV(aPix.y(),aTens.x()*aWeight);
          aSumTy.DIm().AddV(aPix.y(),aTens.y()*aWeight);
     }

     ExpFilterOfStdDev(aPop.DIm()  ,5,mSigmaTensZ1/mZoomRed);
     ExpFilterOfStdDev(aSumTx.DIm(),5,mSigmaTensZ1/mZoomRed);
     ExpFilterOfStdDev(aSumTy.DIm(),5,mSigmaTensZ1/mZoomRed);

     DivImageInPlace(aSumTx.DIm(),aSumTx.DIm(),aPop.DIm());
     DivImageInPlace(aSumTy.DIm(),aSumTy.DIm(),aPop.DIm());


     mImTgt = cIm1D<tREAL8>(mSzRed.y(),nullptr,eModeInitImage::eMIA_Null);
     mDImTgt = & (mImTgt.DIm());

     for(const auto anY : aPop.DIm())
     {
         cPt2dr aTens(aSumTx.DIm().GetV(anY),aSumTy.DIm().GetV(anY));
         cPt2dr aRhoTeta = ToPolar(aTens,0.0);
         tREAL8 aTeta = aRhoTeta.y();
         tREAL8 aTgt = tan(aTeta/2.0);

         mDImTgt->SetV(anY,aTgt);
     }

     mYC = DetectSymByTgt();


     mImIntegr  = cIm1D<tREAL8>(mSzRed.y(),nullptr,eModeInitImage::eMIA_Null);
     mDImIntegr = &(mImIntegr.DIm());

     DoIntegrale(-1,0);
     DoIntegrale(1,mSzRed.y());

}

std::string cAppliNewFrange::NameVisu(const std::string & aPref) const
{
    return mPhProj.DirVisuAppli() + aPref + mNameIm;
}


void cAppliNewFrange::MakeImageMaxLoc()
{
    mImMax = cIm2D<tU_INT1>(mSzRed,nullptr,eModeInitImage::eMIA_Null);
    mDImMax = & (mImMax.DIm());




    int aNbMaxX=8;
    int aNbMaxY=3;

    for (const auto aPix : mDImRed->Interior(1+std::max(aNbMaxX,aNbMaxY)))
    {
         if (     IsMax(aPix,cPt2di(1,0),aNbMaxX)
              ||  IsMax(aPix,cPt2di(0,1),aNbMaxY)
            )
         {
              mDImMax->SetV(aPix,1);
         }
    }
}


void  cAppliNewFrange::DoOneImage(const std::string & aNameIm)
{
    mNameIm = aNameIm;
    mImZ1 = tIm::FromFile(mNameIm);
    mDImZ1 = & (mImZ1.DIm());
    mImRed = mImZ1.BiCubicDeZoom(mZoomRed);
    mDImRed = &(mImRed.DIm());
    mSzRed = mDImRed->Sz();

    mImRedBlur = mImRed.Dup();
    mDImRedBlur = & (mImRedBlur.DIm());
    ExpFilterOfStdDev(*mDImRedBlur,5,2.0,5.0);


    MakeImageMaxLoc();
    /*


    mImMax = cIm2D<tU_INT1>(mSzRed,nullptr,eModeInitImage::eMIA_Null);
    mDImMax = & (mImMax.DIm());




    int aNbMaxX=8;
    int aNbMaxY=3;

    for (const auto aPix : mDImRed->Interior(1+std::max(aNbMaxX,aNbMaxY)))
    {
         if (     IsMax(aPix,cPt2di(1,0),aNbMaxX)
              ||  IsMax(aPix,cPt2di(0,1),aNbMaxY)
            )
         {
              mDImMax->SetV(aPix,1);
         }
    }*/


    if (mDoSimul)
        MakeImSimul();

    if (mDoVisu)
    {
        tREAL8 aNbX=400;

        mImVisu = cRGBImage(mSzRed + cPt2di(aNbX,0),cRGBImage::White);
        tElIm aVMin,aVMax;
        GetBounds(aVMin,aVMax,*mDImRed);

        for (const auto aPix : *mDImRed)
            mImVisu.SetGrayPix(aPix,mDImRed->GetV(aPix)*(255.0/aVMax));
    }


    MakeImTgt();


    if (mDoVisu)
    {
        for (const auto aPix : * mDImMax)
            if (mDImMax->GetV(aPix) != 0 )
                mImVisu.SetRGBPix(aPix,cRGBImage::Cyan);

        for(const auto anY : mImTgt.DIm())
        {
           cPt2di  aPt(round_ni(mSzRed.x()+aNbX/2.0+mDImTgt->GetV(anY)*10.0),anY.x());
           mImVisu.SetRGBPix(cPt2di(mSzRed.x()+aNbX/2,anY.x()),cRGBImage::Green);
           mImVisu.SetRGBPix(aPt,cRGBImage::Red);

           int aXInt = mDImIntegr->GetV(anY.x());
           mImVisu.SetRGBPix(cPt2di(aXInt+100,anY.x()),cRGBImage::Red);
        }

        mImVisu.DrawLine(cPt2dr(0,mYC),cPt2dr(mSzRed.x()+aNbX,mYC),cRGBImage::Blue,1.0);

        mImVisu.ToFile(NameVisu("ImRed"));
        mDImRedBlur->ToFile(NameVisu("Blured"));
     }
}


};  // NS_FrangesDetect


using  namespace NS_FrangesDetect;

// ============================= Old version ===============

tMMVII_UnikPApli Alloc_NewFrangeDetect(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec)
{
      return tMMVII_UnikPApli(new cAppliNewFrange(aVArgs,aSpec));
}


cSpecMMVII_Appli  TheSpecAppliFranges_2
(
     "ExtractFranges_2",
      Alloc_NewFrangeDetect,
      "New new ... Extraction of Franges, version N (image filter)",
      {eApF::ImProc},
      {eApDT::Image},
      {eApDT::Console},
      __FILE__
);



};
