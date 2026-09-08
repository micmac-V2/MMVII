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
        std::string              mFilterAdd;   // pattern to filter GCP by additional info
        bool                     mVerbose;
        bool                     mAdjustSim;   // try to adjust a similarity?

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
    mVerbose      (false),
    mAdjustSim    (false)
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
                << AOpt2007(mAdjustSim, "AdjSim"  ,"Adjust a similarity",{eTA2007::HDV})
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


    std::vector<cPt3dr> aVPts1;
    std::vector<cPt3dr> aVPts2;
    std::vector<std::string> aVPtsNames;


    for (auto& aGCPGnd1: aSetMes1.MesGCP())
    {
        for (auto& aGCPGnd2: aSetMes2.MesGCP())
        {
            if (aGCPGnd1.mNamePt == aGCPGnd2.mNamePt)
            {
                aVPts1.push_back(aGCPGnd1.mPt);
                aVPts2.push_back(aGCPGnd2.mPt);
                aVPtsNames.push_back(aGCPGnd1.mNamePt);
            }
        }
    }

    if (mAdjustSim && (aVPts1.size()>=3))
    {
        cSimilitud3D<tREAL8>  aSim;
        double                mRes2;
        //estimate 3d similarity
        aSim = aSim.StdGlobEstimate(aVPts1,aVPts2,&mRes2,nullptr,cParamCtrlOpt::Default());
        StdOut() << "Similarity residual = " << mRes2 << "\n";
        StdOut() << "Similarity scale = "        << aSim.Scale() << "\n";
        StdOut() << "Similarity translation = "  << aSim.Tr() << "\n";
        StdOut() << "Similarity rotation Mat = "     << aSim.Rot().Mat();
        StdOut() << "Similarity rotation WKP deg = "<< aSim.Rot().ToWPK() * (180./M_PI) <<"\n\n";
        for (auto &aPt2: aVPts2)
            aPt2 = aSim.Inverse(aPt2);
        StdOut() << "Points diff AFTER similarity\n";
    }
    if (mAdjustSim && (mGndCoord2.size()<3))
    {
        StdOut() << "Not enought common point for similarity estimation\n\n";
    }

    cAvgAndBoundVals<tREAL8>  aVDiff;
    cAvgAndBoundVals<tREAL8>  aVDiff2;

    for (size_t i=0; i<aVPts1.size();++i)
    {
        auto aDiff = aVPts2[i] - aVPts1[i];
        auto aNormDiff = Norm2(aDiff);
        if (mVerbose)
            StdOut() << aVPtsNames[i] << ": "<< aNormDiff << "   " << aDiff <<"\n";
        aVDiff.Add(aNormDiff);
        aVDiff2.Add(Square(aNormDiff));
    }
    if (mVerbose)
        StdOut() << "\n";

    StdOut() << mGndCoord2 << " - " << mPhProj.DPGndPt3D().DirIn()
             << " => "<<aVPts1.size()<<" common points, " << Color::descr << "AVG="<< Color::end << sqrt(aVDiff2.Avg())
             << Color::descr <<" RMSE=" << Color::end << sqrt(aVDiff2.Avg()) << Color::descr << " MAX="<< Color::end << aVDiff.VMax() <<"\n";

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

