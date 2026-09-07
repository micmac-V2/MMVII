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



/* ==================================================== */
/*                                                      */
/*                   cAppli_CompareGCP                  */
/*                                                      */
/* ==================================================== */

class cAppli_CompareGCP : public cMMVII_Appli
{
    public :

        cAppli_CompareGCP(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli &);
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;

        std::vector<cOneHelpSampleCmp> Samples() const override;


    private :
        cPhotogrammetricProject  mPhProj;

        std::string              mGndCoord2;
        std::string              mFilterName;  // pattern to filter names of GCP
        std::string              mFilterAdd;  // pattern to filter GCP by additional info
        bool                     mVerbose;

};

cAppli_CompareGCP::cAppli_CompareGCP
(
    const std::vector<std::string> &  aVArgs,
    const cSpecMMVII_Appli & aSpec
) :
    cMMVII_Appli  (aVArgs,aSpec),
    mPhProj       (*this),
    mFilterName   (),
    mFilterAdd    (),
    mVerbose      (false)
{
}



cCollecSpecArg2007 & cAppli_CompareGCP::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return anArgObl
                << mPhProj.DPGndPt3D().ArgDirInMand()
                << Arg2007(mGndCoord2,"Second ground coord folder",{{eTA2007::Input},{eTA2007::ObjCoordWorld}})
            ;
}

cCollecSpecArg2007 & cAppli_CompareGCP::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    return anArgOpt
               << AOpt2007(mFilterName, "Filter", "Pattern to filter GCP by name")
               << AOpt2007(mFilterAdd, "FilterAdd", "Pattern to filter GCP by additional info")
               << AOpt2007(mVerbose, "Verbose"  ,"Do show each difference",{eTA2007::HDV})
            ;

}

std::vector<cOneHelpSampleCmp>  cAppli_CompareGCP::Samples() const
{
    return {
            {"MMVII ReportGCPCmp Init Final Filter='.*_'"}
    };
}




int cAppli_CompareGCP::Exe()
{
    mPhProj.FinishInit();

    cDirsPhProj &  aDirGndCoord2 = *mPhProj.NewDPIn(eTA2007::ObjCoordWorld,mGndCoord2);

    MMVII_INTERNAL_ASSERT_tiny(IsDirectory(mPhProj.DPGndPt3D().FullDirIn()),
                               "CoordWolrd " + mPhProj.DPGndPt3D().FullDirIn() + " does not exist");

    MMVII_INTERNAL_ASSERT_tiny(IsDirectory(aDirGndCoord2.FullDirIn()),
                               "CoordWolrd " + mGndCoord2 + " does not exist");

    cSetMesGndPt  aSetMes1;
    mPhProj.LoadGCP3D(aSetMes1,nullptr,"",mFilterName,mFilterAdd);

    cSetMesGndPt  aSetMes2;
    mPhProj.LoadGCP3DFromFolder(aDirGndCoord2.FullDirIn(), aSetMes2, nullptr, "", mFilterName, mFilterAdd);


    cAvgAndBoundVals<tREAL8>  aVDiff;
    for (auto& aGCPGnd1: aSetMes1.MesGCP())
    {
        for (auto& aGCPGnd2: aSetMes2.MesGCP())
        {
            if (aGCPGnd1.mNamePt == aGCPGnd2.mNamePt)
            {
                auto aDiff = aGCPGnd2.mPt - aGCPGnd1.mPt;
                auto aNormDiff = Norm2(aDiff);
                if (mVerbose)
                    StdOut() << aGCPGnd1.mNamePt << ": "<< aNormDiff << "   " << aDiff <<"\n";
                aVDiff.Add(aNormDiff);
            }
        }
    }

    StdOut() << "\n"<< mGndCoord2 << " - " << mPhProj.DPGndPt3D().DirIn()
             <<" =>   Average: "<< aVDiff.Avg() << "   Max: "<< aVDiff.VMax() <<"\n";

   return EXIT_SUCCESS;
}

/* ==================================================== */
/*                                                      */
/*               MMVII                                  */
/*                                                      */
/* ==================================================== */



tMMVII_UnikPApli Alloc_CompareGCP(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppli_CompareGCP(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_CompareGCP
(
     "ReportGCPCmp",
      Alloc_CompareGCP,
      "Compare GCP 3D coords",
      {eApF::GCP},
      {eApDT::ObjCoordWorld},
      {eApDT::Xml},
      __FILE__
);



}; // MMVII

