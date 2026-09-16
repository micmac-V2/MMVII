#include "MMVII_Ptxd.h"
#include "MMVII_SysSurR.h"
#include "MMVII_Sensor.h"
#include "MMVII_PCSens.h"
#include "MMVII_Tpl_Images.h"
#include "MMVII_BundleAdj.h"


namespace MMVII
{


template<class Type>  cIm2D<tU_INT1>  Image8Bits (const cDataIm2D<Type> & anIm)
{
     cIm2D<tU_INT1> aRes = anIm.Sz();
     Type aVMin,aVMax;

    GetBounds(aVMin,aVMax,anIm);

    MMVII_INTERNAL_ASSERT_always(aVMax!=0,"Image8Bits Vmax=0");

    for (int aK=1 ; aK<anIm.NbElem() ; aK++)
    {
        const Type & aV = anIm.GetRDL(aK);
        aRes.DIm().GetRDL(aK)  = aV * (255.0 / aVMax);

    }

     return aRes;
}


/* ==================================================== */
/*                                                      */
/*              cAppli_StackIm                          */
/*                                                      */
/* ==================================================== */

class cAppli_StackIm : public cMMVII_Appli
{
     public :

        typedef cIm2D<tREAL4> tIm;
        typedef cDataIm2D<tREAL4> tDIm;

        cAppli_StackIm(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli &);
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;

     private :

        std::string  mSpecImIn;
        std::string  mNameOut;
        bool          m8Bits;
        tREAL8        mFreq;
};

cAppli_StackIm::cAppli_StackIm(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec):
        cMMVII_Appli   (aVArgs,aSpec),
        mNameOut       ("Stack.tif"),
        m8Bits         (false),
        mFreq          (1e-1)
{
}

cCollecSpecArg2007 & cAppli_StackIm::ArgObl(cCollecSpecArg2007 & anArgObl)
{
      return anArgObl
              << Arg2007(mSpecImIn,"Pattern/file for images",{{eTA2007::MPatFile,"0"},{eTA2007::FileDirProj}})


/*
              <<  mPhProj.DPPointsMeasures().ArgDirInMand()
              <<  mPhProj.DPOrient().ArgDirOutMand()
*/
           ;
}

cCollecSpecArg2007 & cAppli_StackIm::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{

    return anArgOpt

            << AOpt2007(mFreq,"FreqShow","Frequence where we show how many remain",{eTA2007::HDV})
            << AOpt2007(mNameOut,CurOP_Out,"Destimation file",{eTA2007::HDV})
            << AOpt2007(m8Bits,"Do8Bits","Generate 8 bits image",{eTA2007::HDV})

           ;
}


int cAppli_StackIm::Exe()
{
     //for (const auto &  aNameIm : VectMainSet(0))

    auto aVecIm = VectMainSet(0);

    tIm aIm0 = tIm::FromFile(aVecIm.at(0));
    tDIm & aDIm0 = aIm0.DIm();
    aDIm0.InitCste(0);

    for (size_t aKIm=0 ; aKIm<aVecIm.size() ; aKIm++)
    {
         const auto & aNameIm  = aVecIm.at(aKIm);
         tIm aImK = tIm::FromFile(aNameIm);

         AddIn(aDIm0,aImK.DIm());

         if (SignalAtFrequence(aKIm,mFreq,0.0))
             StdOut() << " StackIm still to do " << aVecIm.size() - aKIm << "\n";
    }
    DivCsteIn(aDIm0,tREAL4(aVecIm.size()));

    if (m8Bits)
        Image8Bits(aDIm0).DIm().ToFile(mNameOut);
    else
        aDIm0.ToFile(mNameOut);


    return EXIT_SUCCESS;
};




/* ==================================================== */
/*                                                      */
/*                                                      */
/*                                                      */
/* ==================================================== */


tMMVII_UnikPApli Alloc_StackIm(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppli_StackIm(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_StackIm
(
     "ImageStack",
      Alloc_StackIm,
      "Stack a serie of images",
      {eApF::ImProc},
      {eApDT::Image},
      {eApDT::Image},
      __FILE__
);



}; // MMVII




