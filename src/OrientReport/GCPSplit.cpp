#include "MMVII_DeclareAllCmd.h"
#include "MMVII_Sensor.h"

namespace MMVII
{

/* ********************************************************** */
/*                                                            */
/*                 cAppli_GCPSplit                           */
/*                                                            */
/* ********************************************************** */
class cAppli_GCPSplit : public cMMVII_Appli
{
    public:
        cAppli_GCPSplit(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec);
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;

    private:

        cPhotogrammetricProject   mPhProj;
        double                    mSplitRatio;
        std::string               mOutPrefix;
        bool                      mVerbose;

};

cAppli_GCPSplit::cAppli_GCPSplit(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
    cMMVII_Appli (aVArgs,aSpec),
    mPhProj      (*this),
    mSplitRatio  (0.5),
    mOutPrefix   (""),
    mVerbose(false)
{}

cCollecSpecArg2007 & cAppli_GCPSplit::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return anArgObl
           << mPhProj.DPGndPt3D().ArgDirInOpt()
           << Arg2007(mSplitRatio,"Split ratio",{eTA2007::HDV});

}

cCollecSpecArg2007 & cAppli_GCPSplit::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    return anArgOpt
        //   <<  mPhProj.DPGndPt3D().ArgDirOutOpt()
           <<  AOpt2007(mOutPrefix,"OutPrefix","Output name prefix",{eTA2007::HDV})
           <<  AOpt2007(mVerbose,"Verbose","Verbose print",{eTA2007::HDV,eTA2007::Tuning});
}

int cAppli_GCPSplit::Exe()
{
    mPhProj.FinishInit();

    cSetMesGndPt aSetMes;
    cSetMesGnd3D aSetMesGCP("GCPSplit");
    cSetMesGnd3D aSetMesCP("GCPSplit");

    mPhProj.LoadGCP3D(aSetMes);

    std::vector<cMes1Gnd3D>   aSetMesVec = aSetMes.AllMesGCP().Measures();
    size_t                    aNbPts = aSetMesVec.size();
    size_t                    aNbPtsGCP = int(mSplitRatio*aNbPts);


    std::vector<double> aMinDistVec(aNbPts,std::numeric_limits<double>::max());
    std::vector<int> aSelPts;


    // add first randomly selected points
    int aFirstRandPt = RandUnif_N(aNbPts);
    int aCur = aFirstRandPt;
    aSetMesGCP.AddMeasure3D(cMes1Gnd3D(aSetMesVec.at(aCur)));


    // iterate until you reach the expected num of points
    while (aSetMesGCP.Measures().size() < aNbPtsGCP)
    {
        // update min distances to the current point
        int    aMaxDistId=1e5;
        double aMaxDist =0;
        cPt3dr aCurPt = aSetMesVec.at(aCur).mPt;
        for (size_t aK=0; aK<aNbPts; aK++)
        {
            aMinDistVec.at(aK) = std::min(aMinDistVec.at(aK), SqN2( aCurPt-aSetMesVec.at(aK).mPt ));
            if (aMinDistVec.at(aK) > aMaxDist)
            {
                aMaxDistId = aK;
                aMaxDist = aMinDistVec.at(aK);
            }
        }

        // update measurement structure with the farthest point
        cMes1Gnd3D aSelGCP (aSetMesVec.at(aMaxDistId));
        aSetMesGCP.AddMeasure3D(aSelGCP);

        if (mVerbose)
            StdOut() << "GCP " << aSelGCP.mNamePt << " " << aSelGCP.mPt << std::endl;


        // update the index of current points
        aCur = aMaxDistId;

    }

    // collect check points (non-gcp points)
    int aK=0;
    for (auto aDist : aMinDistVec)
    {
        if (aDist!=0 && aK!=aFirstRandPt)
        {
            aSetMesCP.AddMeasure3D(cMes1Gnd3D(aSetMesVec.at(aK)));
            if (mVerbose)
                StdOut() << "CP " << aSetMesVec.at(aK).mNamePt << " " << aSetMesVec.at(aK).mPt << std::endl;
        }
        aK++;
    }

    std::string aNameGCP = "GCP_" + mOutPrefix + mPhProj.DPGndPt3D().DirIn();
    std::string aNameCP =  "CP_" + mOutPrefix + mPhProj.DPGndPt3D().DirIn();


    mPhProj.SaveGCP3D(aSetMesGCP,aNameGCP,true);
    mPhProj.SaveGCP3D(aSetMesCP,aNameCP,true);

    StdOut() << " * GCP " << aSetMesGCP.Measures().size() << "/" << aSetMes.AllMesGCP().Measures().size() << ": " << aNameGCP << std::endl;
    StdOut() << " * CP " << aSetMesCP.Measures().size() << "/" << aSetMes.AllMesGCP().Measures().size() << ": " << aNameCP << std::endl;

    return EXIT_SUCCESS;
}


/* ====================================================== */
/*                     GCPSplit                           */
/* ====================================================== */

tMMVII_UnikPApli Alloc_GCPSplit(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
    return tMMVII_UnikPApli(new cAppli_GCPSplit(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_GCPSplit
    (
        "GCPSplit",
        Alloc_GCPSplit,
        "Split a set of GCPs into check and control pts",
        {eApF::GCP},
        {eApDT::ObjCoordWorld},
        {},
        __FILE__
    );


}; // MMVII
