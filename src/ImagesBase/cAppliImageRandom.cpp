#include "MMVII_Tpl_Images.h"
#include "MMVII_Linear2DFiltering.h"
#include "MMVII_DeclareCste.h"
#include "MMVII_Geom2D.h"
#include "MMVII_Interpolators.h"
#include "MMVII_Tpl_ElemStrToVal.h"


namespace MMVII
{

struct cParamNoise
{
    tREAL8 mSigma ;
    tREAL8 mWeight  ;


    ARG2007_STRUCT_FIELDS (
        mWeight, FieldSem({{eTA2007::AddCom,"Weight on noise"}}),
        mSigma,FieldSem({{eTA2007::AddCom,"Sigma on noise"}})
    )

};




struct cParamComb
{
    tREAL8 mPeriod ;
    tREAL8 mSigma ;
    tREAL8 mIntervAmpl ;



    ARG2007_STRUCT_FIELDS (
        mPeriod,FieldSem({{eTA2007::AddCom,"Period of diracs"}}),
        mSigma, FieldSem({{eTA2007::AddCom,"Sigma applied to dirac"}}),
        mIntervAmpl, FieldSem({{eTA2007::AddCom,"Random interval for each dirac"}})
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


       tIm       mIm;
       tDIm *    mDIm;
        
};  





cAppliGenRandomImage::cAppliGenRandomImage
(
      const std::vector<std::string> & aVArgs,
      const cSpecMMVII_Appli & aSpec
) :
   cMMVII_Appli (aVArgs,aSpec) ,
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
    return
          anArgOpt

      // << cHeaderSectionArg("Topo")

      << AOpt2007 ( mParamNoise, "GaussNoise", "Parameter for gaussian noise",{eTA2007::CanRepeat})
      << AOpt2007 ( mParamCombs, "Combs", "Parameter for combs",{eTA2007::CanRepeat})


   ;
}



int cAppliGenRandomImage::Exe()
{
    mIm = tIm(mSz,nullptr,eModeInitImage::eMIA_Null);
    mDIm = & mIm.DIm();
    tREAL8 aSumW = 0.0;


    for (const auto & aGN : mParamNoise)
    {
       tIm aImRand(mSz,nullptr,eModeInitImage::eMIA_RandCenter);
       ExpFilterOfStdDev(aImRand.DIm(),5,aGN.mSigma);

       NormalizedAvgDev(aImRand.DIm(),1e-5);

       AddMulImageCsteInPlace(*mDIm,aImRand.DIm(),aGN.mWeight);

       aSumW += aGN.mWeight;
       StdOut() <<  " Done Gaussian noise" << aGN.mSigma << " " << aGN.mWeight << "\n";
    }



    for (const auto & aParamC : mParamCombs)
    {
       tIm aImComb(mSz,nullptr,eModeInitImage::eMIA_Null);

       cPt2dr aPer (aParamC.mPeriod,aParamC.mPeriod);
       cPt2dr aIntR  (-aParamC.mIntervAmpl,aParamC.mIntervAmpl);
       cPt2dr aSigma(aParamC.mSigma,aParamC.mSigma);
       cPt2dr aPt;
       for (aPt.x()=aPer.x()/2.0 ;  aPt.x()<mSz.x() ; aPt.x()+=aPer.x())
       {
           for (aPt.y()=aPer.y()/2.0 ;  aPt.y()<mSz.y() ; aPt.y()+=aPer.y())
           {
               if (aImComb.DIm().InsideBL(aPt))
                   aImComb.DIm().AddVBL(aPt,RandInInterval(aIntR.x(),aIntR.y()));
           }
       }
       ExpFilterOfStdDev(aImComb.DIm(),5,aSigma.x(),aSigma.y());

       AddIn(*mDIm,aImComb.DIm());

    //   ExpFilterOfStdDev(aImRand.DIm(),5,aGN.mSigma);

       /*NormalizedAvgDev(aImRand.DIm(),1e-5);

       AddMulImageCsteInPlace(*mDIm,aImRand.DIm(),aGN.mWeight);

       aSumW += aGN.mWeight;
       StdOut() <<  " Done Gaussian noise" << aGN.mSigma << " " << aGN.mWeight << "\n";*/
    }




    mDIm->ToFile("toto.tif");

    StdOut() << " SUMW " << aSumW << "\n";

    return EXIT_SUCCESS;
}


tMMVII_UnikPApli Alloc_cAppliGenRandomImage(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppliGenRandomImage(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_cAppliGenRandomImage
(
     "ImageGenRandom",
      Alloc_cAppliGenRandomImage,
      "Bundle adjusment between images, using several observations/constraint",
      {eApF::Ori},
      {eApDT::Orient},
      {eApDT::Orient},
      __FILE__
);




}; // namespace MMVII

