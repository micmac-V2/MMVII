#include "MMVII_Ptxd.h"
#include "cMMVII_Appli.h"
#include "MMVII_Geom3D.h"
#include "MMVII_PCSens.h"
#include "MMVII_Tpl_Images.h"
#include "MMVII_StaticLidar.h"

/**
   \file GCPCompare.cpp


 */

namespace MMVII
{

/*cSet2D3D

cPerspCamIntrCalib*/
//cSet2D3D

/* ==================================================== */
/*                                                      */
/*                   cAppli_CmpCalib                    */
/*                                                      */
/* ==================================================== */

class cAppli_CmpCalib : public cMMVII_Appli
{
    public :

        cAppli_CmpCalib(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli &,bool isModeLocal);
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;

        std::vector<cOneHelpSampleCmp> Samples() const override;        

    private :

        void CmpCalib(const cPerspCamIntrCalib & aCal1,const cPerspCamIntrCalib & aCal2);

        cPhotogrammetricProject  mPhProj;
        bool                     mIsModeLocal;
        std::string              mNameIm1;
        std::string              mNameIm2;

        std::string              mOri1;
        std::string              mOri2;
};

cAppli_CmpCalib::cAppli_CmpCalib
(
    const std::vector<std::string> &  aVArgs,
    const cSpecMMVII_Appli & aSpec,
    bool  isModeLocal
) :
    cMMVII_Appli  (aVArgs,aSpec),
    mPhProj       (*this),
    mIsModeLocal  (isModeLocal),
    mNameIm1      (MMVII_NONE),
    mNameIm2      (MMVII_NONE),
    mOri1         (MMVII_NONE),
    mOri2         (MMVII_NONE)
{
}



cCollecSpecArg2007 & cAppli_CmpCalib::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    if (mIsModeLocal)
    {
        anArgObl
                << Arg2007(mNameIm1,"Name of first image",{{eTA2007::FileImage},{eTA2007::FileDirProj}})
                << mPhProj.DPOrient().ArgDirInMand()
        ;
    }
    else
    {

    }

    return anArgObl;
}

cCollecSpecArg2007 & cAppli_CmpCalib::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    if (mIsModeLocal)
    {
        anArgOpt
              << AOpt2007(mOri2,"Ori2","Second orientation folder when != first",{{eTA2007::Input},{eTA2007::Orient}})
              << AOpt2007(mNameIm2,"Im2","Second image when != first ",{{eTA2007::Input},{eTA2007::FileImage}})
            ;
     }
    else
    {
    }

     return anArgOpt;
}

std::vector<cOneHelpSampleCmp>  cAppli_CmpCalib::Samples() const
{
    return {
            // {"MMVII ReportGCPCmp Init Final Filter='.*_'"}
    };
}


void cAppli_CmpCalib::CmpCalib(const cPerspCamIntrCalib & aCal1,const cPerspCamIntrCalib & aCal2)
{
    StdOut()  << " O1=" << mOri1 << " F1=" << aCal1.F()
              << " O2=" << mOri2 << " F2=" << aCal2.F()<< "\n";
}


int cAppli_CmpCalib::Exe()
{
   mPhProj.FinishInit();


  if (mIsModeLocal)
  {
      mOri1 = mPhProj.DPOrient().DirIn();
      if (!IsInit(&mOri2))
          mOri2 = mOri1;

      if (!IsInit(&mNameIm2))
          mNameIm2 = mNameIm1;

      cPerspCamIntrCalib * aCal1 = mPhProj.InternalCalibFromStdName(mNameIm1);
      cPerspCamIntrCalib * aCal2 = mPhProj.InternalCalibFromFolderStdName(mOri2,mNameIm2);


      CmpCalib(*aCal1,*aCal2);


  }
  else
  {
  }

   return EXIT_SUCCESS;
}

/* ==================================================== */
/*                                                      */
/*               MMVII                                  */
/*                                                      */
/* ==================================================== */



tMMVII_UnikPApli Alloc_CmpCalib_Local(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppli_CmpCalib(aVArgs,aSpec,true));
}

cSpecMMVII_Appli  TheSpec_CmpCalib_Local
(
     "OriCmpCalibLocal",
      Alloc_CmpCalib_Local,
      "Compare internal calibrations, coming from local folder",
      {eApF::GCP},
      {eApDT::ObjCoordWorld},
      {eApDT::Xml},
      __FILE__
);



}; // MMVII

