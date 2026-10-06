#include "MMVII_Tpl_Images.h"
#include "MMVII_Linear2DFiltering.h"
#include "MMVII_DeclareCste.h"
#include "MMVII_Geom2D.h"
#include "MMVII_Interpolators.h"
#include "MMVII_Tpl_ElemStrToVal.h"


namespace MMVII
{

static const   std::string TheNameCmd("ImageGenRandom");


struct cParamNoise
{
    tREAL8 mSigma ;
    tREAL8 mWeight  ;


    ARG2007_STRUCT_FIELDS (
        mSigma,FieldSem({{eTA2007::AddCom,"Sigma on noise"}}),
        mWeight, FieldSem({{eTA2007::AddCom,"Weight on noise"}})
    )

};




struct cParamComb // peigne ...
{
    tREAL8 mPeriod ;
    tREAL8 mSigma ;
    tREAL8 mWeight ;



    ARG2007_STRUCT_FIELDS (
        mPeriod,FieldSem({{eTA2007::AddCom,"Period of diracs"}}),
        mSigma, FieldSem({{eTA2007::AddCom,"Sigma applied to dirac"}}),
        mWeight, FieldSem({{eTA2007::AddCom,"Ampl of global"}})
    )

};

class cAppliGenRandomImage : public cMMVII_Appli
{
     public :
        cAppliGenRandomImage(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec);
     private :
        //================= typedef part ===============================


        typedef tREAL8         tEl;
        typedef cIm2D<tEl >    tIm;
        typedef cDataIm2D<tEl> tDIm;

        //================================================================
        //       METHODS DECLARATION
        //================================================================

            //------------------------- overidding cMMVII_Appli ----------------
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;
        // std::vector<cOneHelpSampleCmp>  Samples() const override;
        //virtual ~cAppliGenRandomImage();

     private :

       cPt2di                     mSz;
       std::vector<cParamNoise>   mParamNoise;
       std::vector<cParamComb>    mParamCombs;
       std::string                mNameOut;
       eTyNums                    mTypeOut;
       tREAL8                     mOffset;

       tIm       mIm;
       tDIm *    mDIm;
        
};  





cAppliGenRandomImage::cAppliGenRandomImage
(
      const std::vector<std::string> & aVArgs,
      const cSpecMMVII_Appli & aSpec
) :
   cMMVII_Appli (aVArgs,aSpec) ,
   mNameOut     (TheNameCmd + std::string(".tif")),
   mTypeOut     (eTyNums::eTN_REAL4),
   mOffset      (0.0),
   mIm          (cPt2di(1,1)),
   mDIm         (nullptr)
{
}


cCollecSpecArg2007 & cAppliGenRandomImage::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return anArgObl
              << Arg2007(mSz,"Size of images")
           ;
}



cCollecSpecArg2007 & cAppliGenRandomImage::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    return anArgOpt
            << AOpt2007 ( mParamNoise, "GaussNoise", "Parameter for gaussian noise",{eTA2007::CanRepeat})
            << AOpt2007 ( mParamCombs, "Combs", "Parameter for combs (aka \"peignes\")",{eTA2007::CanRepeat})
            << AOpt2007 ( mNameOut, CurOP_Out , "Output file",{eTA2007::HDV})
            << AOpt2007 ( mTypeOut, "Type" , "Type of element of  generated file",{eTA2007::HDV})
            << AOpt2007 ( mOffset,"Offset","Offset to add, def depend of type")
   ;
}



int cAppliGenRandomImage::Exe()
{
    mIm = tIm(mSz,nullptr,eModeInitImage::eMIA_Null);
    mDIm = & mIm.DIm();

    if (! IsInit(&mNameOut))
    {
        mNameOut = mSpecs.Name() + ".tif";
    }

    const cVirtualTypeNum & aVTN = cVirtualTypeNum::FromEnum(mTypeOut);
    if (! IsInit(&mOffset))
    {
        mOffset = aVTN.CenteredValue();
    }
    mDIm->InitCste(mOffset);


    for (const auto & aGN : mParamNoise)
    {
       tIm aImRand(mSz,nullptr,eModeInitImage::eMIA_RandCenter);
       ExpFilterOfStdDev(aImRand.DIm(),5,aGN.mSigma);

       GenNormalizedAvgDev(aImRand.DIm(),1e-20);

       AddMulImageCsteInPlace(*mDIm,aImRand.DIm(),aGN.mWeight);

       StdOut() <<  " Done Gaussian noise" << aGN.mSigma << " " << aGN.mWeight << "\n";
    }



    for (const auto & aParamC : mParamCombs)
    {
       tIm aImComb(mSz,nullptr,eModeInitImage::eMIA_Null);

       cPt2dr aPer (aParamC.mPeriod,aParamC.mPeriod);
       cPt2dr aSigma(aParamC.mSigma,aParamC.mSigma);
       cPt2dr aPt;

       // make sum of dirac with random ampl
       for (aPt.x()=aPer.x()/2.0 ;  aPt.x()<mSz.x() ; aPt.x()+=aPer.x())
       {
           for (aPt.y()=aPer.y()/2.0 ;  aPt.y()<mSz.y() ; aPt.y()+=aPer.y())
           {
               if (aImComb.DIm().InsideBL(aPt))
               {
                   aImComb.DIm().AddVBL(aPt,RandInInterval(-1,1));
               }
           }
       }
       // convoluate
       ExpFilterOfStdDev(aImComb.DIm(),5,aSigma.x(),aSigma.y());
       // normalize
       GenNormalizedAvgDev(aImComb.DIm(),1e-20,1.0);
       // add
       AddMulImageCsteInPlace(*mDIm,aImComb.DIm(),aParamC.mWeight);
    }

    mDIm->ToFile(mNameOut,mTypeOut);


    return EXIT_SUCCESS;
}


tMMVII_UnikPApli Alloc_cAppliGenRandomImage(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppliGenRandomImage(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_cAppliGenRandomImage
(
      TheNameCmd,
      Alloc_cAppliGenRandomImage,
      "Bundle adjusment between images, using several observations/constraint",
      {eApF::Ori},
      {eApDT::Orient},
      {eApDT::Orient},
      __FILE__
);




}; // namespace MMVII

