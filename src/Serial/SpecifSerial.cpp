#include "MMVII_Bench.h"
#include "MMVII_Class4Bench.h"
#include "MMVII_2Include_Serial_Tpl.h"

#include "MMVII_Geom2D.h"
#include "MMVII_PCSens.h"
#include "Serial.h"
#include "MMVII_MeasuresIm.h"


/** \file DemoSerial.cpp
    \brief file for generating spec+samples of serialization

*/

namespace MMVII
{

/* *********************************************************** */
/*                                                             */
/*                   cAppliSpecSerial                          */
/*                                                             */
/* *********************************************************** */

class cAppliSpecSerial : public cMMVII_Appli
{
     public :

        cAppliSpecSerial(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec);

     private :
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;
        std::string mOutDir;

};


cAppliSpecSerial::cAppliSpecSerial
(
    const std::vector<std::string> & aVArgs,
    const cSpecMMVII_Appli & aSpec
) :
    cMMVII_Appli   (aVArgs,aSpec),
    mOutDir(DirRessourcesMMVII() + "SpecifSerial/")

{
}

cCollecSpecArg2007 & cAppliSpecSerial::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return anArgObl
   ;
}


cCollecSpecArg2007 & cAppliSpecSerial::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
   return anArgOpt
           << AOpt2007(mOutDir,"OutDir","Directory for generated specifiction files", {eTA2007::HDV})
          ;
}


void GenSpec_BitEncoding(const std::string & aDir);
void GenSpec_SysCo(const std::string & aDir);

template<typename T>
static inline void GenSpec(const std::string & aFile)
{
    SpecificationSaveInFile<T>(aFile +".xml");
    SpecificationSaveInFile<T>(aFile +".json");

}

int  cAppliSpecSerial::Exe()
{
   CreateDirectories(mOutDir,false);
   mOutDir += "/";
   GenSpec_BitEncoding(mOutDir);
   GenSpec_SysCo(mOutDir);
   GenSpec<cSysCoData>(mOutDir+"SysCo");
   GenSpec<tNameSet>(mOutDir+"SetName");
   GenSpec<cSetMesPtOf1Im>(mOutDir+"SetMesureIm");
   GenSpec<cSetMesGnd3D>(mOutDir+"SetMesureGCP");
   GenSpec<cComputeAssociation>(mOutDir+"ComputeAssociation");
   GenSpec<cCamDataBase>(mOutDir+"CamDataBase");


   return EXIT_SUCCESS;
}

/* *********************************************************** */
/*                                                             */
/*                           ::                                */
/*                                                             */
/* *********************************************************** */


tMMVII_UnikPApli Alloc_cAppliSpecSerial(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppliSpecSerial(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_SpecSerial
(
     "GenerateSpecifSerial",
      Alloc_cAppliSpecSerial,
      "Generate specification+some sample for serialization",
      {eApF::Project,eApF::Test},
      {eApDT::None},
      {eApDT::Xml},
      __FILE__
);

/*
*/




};
