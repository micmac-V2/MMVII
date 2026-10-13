//#include "MMVII_PCSens.h"

#include "MMVII_ExtractLines.h"
#include "MMVII_Sensor.h"

#include "cAppliIA.h"


namespace MMVII
{

/* =============================================== */
/*                                                 */
/*                 cAppliBubbles                   */
/*                                                 */
/* =============================================== */


/**  An application for  computing line extraction. Created for CERN
 *  porject, and parametrization is more or less optimized for this purpose.
 */

//    cIm2D<tIm> anIm = cIm2D<tIm>::FromFile(mNameCurIm);

class cAppli_IA : public cMMVII_Appli
{
     public :
        cAppli_IA(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec);
     private :
        typedef tREAL4            tElIm;
        typedef cIm2D<tElIm >     tIm;
        typedef cDataIm2D<tElIm > tDIm;

        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;
         std::vector<cOneHelpSampleCmp>  Samples() const override;
        virtual ~cAppli_IA();

         std::string mFileTaux;



};


cAppli_IA::cAppli_IA(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
    cMMVII_Appli      (aVArgs,aSpec)
{

}

cAppli_IA::~cAppli_IA()
{
}

cCollecSpecArg2007 & cAppli_IA::ArgObl(cCollecSpecArg2007 & anArgObl)
{
      return    anArgObl
           // <<  Arg2007(mPatImage,"Name of input Image", {eTA2007::FileDirProj,{eTA2007::MPatFile,"0"}})
      ;
}

cCollecSpecArg2007 & cAppli_IA::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    return anArgOpt
             << AOpt2007(mFileTaux,"FileTaux","File for change conversion",{eTA2007::HDV})
             /*
            << AOpt2007(mDerFact,"DerFact","Deriche factor for graident",{eTA2007::HDV})*/


            ;
}

std::vector<cOneHelpSampleCmp>  cAppli_IA::Samples() const
{
   return {
        };
}








int cAppli_IA::Exe()
{

    CONVERTIR(mFileTaux);


    return EXIT_SUCCESS;
}


tMMVII_UnikPApli Alloc_AppiIA(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec)
{
      return tMMVII_UnikPApli(new cAppli_IA(aVArgs,aSpec));
}


cSpecMMVII_Appli  TheSpecAppliIA
(
     "IA",
      Alloc_AppiIA,
      "Extraction of Bubbles",
      {eApF::ImProc},
      {eApDT::Image},
      {eApDT::Console},
      __FILE__
);

#if (0)
#endif

};
