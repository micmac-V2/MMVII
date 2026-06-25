#include "MMVII_DeclareAllCmd.h"
#include "MMVII_Sensor.h"

namespace MMVII
{


/* ********************************************************** */
/*                                                            */
/*                 cAppli_OriRel2Im                           */
/*                                                            */
/* ********************************************************** */

class cAppli_OriRelPipeline : public cMMVII_Appli
{
    public :

        cAppli_OriRelPipeline(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec);
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;

    private :

        cPhotogrammetricProject   mPhProj;
        std::string               mPatIm;
     //   std::string               mOriRelOut;
     //   std::string               mOriGlobOut;
        bool                      mExe;
};

cAppli_OriRelPipeline::cAppli_OriRelPipeline(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
    cMMVII_Appli  (aVArgs,aSpec),
    mPhProj (*this),
//    mOriRelOut("Relative"),
//    mOriGlobOut("Global"),
    mExe(true)
{
}

cCollecSpecArg2007 & cAppli_OriRelPipeline::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return anArgObl
           << Arg2007(mPatIm,"Pattern of images",{{eTA2007::MPatFile,"0"},{eTA2007::FileDirProj}})
           << mPhProj.DPOriRel().ArgDirOutMand()
           << mPhProj.DPOrient().ArgDirInMand("Input calibration");

}

cCollecSpecArg2007 & cAppli_OriRelPipeline::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    return anArgOpt
           <<  mPhProj.DPOrient().ArgDirOutOpt()
           <<  mPhProj.DPMulTieP().ArgDirInOpt()
           <<  mPhProj.DPTieP().ArgDirInOpt()
           <<  mPhProj.DPGndPt2D().ArgDirInOpt()
           <<  mPhProj.DPMulTieP().ArgDirOutOpt("VirTP","Output folder for virtual tie points")
           <<  AOpt2007(mExe,"Exe","Execute pipeline",{eTA2007::HDV});
}

int cAppli_OriRelPipeline::Exe()
{
    mPhProj.FinishInit();

    MMVII_INTERNAL_ASSERT_always((mPhProj.DPTieP().DirInIsInit() ||
                                  mPhProj.DPMulTieP().DirInIsInit() ||
                                  mPhProj.DPGndPt2D().DirInIsInit()),"Input tie-points not provided");



    std::list<cParamCallSys>  aLPCS;

    // ===============================================
    // select pairs
    std::string aTieP = mPhProj.DPTieP().DirInIsInit() ? ("InTieP="+mPhProj.DPTieP().DirIn()) :
                    (mPhProj.DPMulTieP().DirInIsInit() ? ("InMulTieP="+mPhProj.DPMulTieP().DirIn()) :
                                                         ("InObjMesInstr="+mPhProj.DPGndPt2D().DirIn()) );
    cParamCallSys aComSelPairs(cMMVII_Appli::FullBin(),
                               "OriPoseSelecAllPAir",
                               mPatIm,
                               mPhProj.DPOriRel().DirOut(),
                               aTieP);
    aLPCS.push_back(aComSelPairs);



    // ===============================================
    // compute relative orientations between pairs
    cParamCallSys aComEstimRelPair(cMMVII_Appli::FullBin(),
                                   "OriPoseEstimRelAllPairs");

    for (size_t aKP=3; aKP<mArgv.size() ; aKP++)
    {
        if (starts_with(mArgv[aKP], "Exe"))
            continue;
        if (starts_with(mArgv[aKP], "VirTP")) //omit virtual points
            continue;
        if (starts_with(mArgv[aKP], "OutOri")) // omit output ori
            continue;

        aComEstimRelPair.AddArgs(mArgv[aKP]);
    }
    aLPCS.push_back(aComEstimRelPair);



    // ===============================================
    // compute relative orientations within triplets
    cParamCallSys aComEstimRelTri(cMMVII_Appli::FullBin(),
                                  "OriPoseEstimRelAllTriplets");


    for (size_t aKP=3; aKP<mArgv.size() ; aKP++)
    {
        if (starts_with(mArgv[aKP], "Exe"))
            continue;
        if (starts_with(mArgv[aKP], "OutOri")) // omit output ori
            continue;

        aComEstimRelTri.AddArgs(mArgv[aKP]);
    }
    aLPCS.push_back(aComEstimRelTri);


    // ===============================================
    // global hierarchical SFM
    if (mPhProj.DPOrient().DirOutIsInit())
    {
        cParamCallSys aComHSfM (cMMVII_Appli::FullBin(),
                               "OriPoseEstimGlobHierarchical",
                               mPatIm);

        // add relative and calibration
        for (size_t aKP=3; aKP<5 ; aKP++)
        {
            aComHSfM.AddArgs(mArgv[aKP]);
        }

        // add output global ori
        aComHSfM.AddArgs(mPhProj.DPOrient().DirOut());

        for (size_t aKP=5; aKP<mArgv.size() ; aKP++)
        {
            if (starts_with(mArgv[aKP], "InMulTieP"))
            {
                // use virtual points
                if (mPhProj.DPMulTieP().DirOutIsInit())
                {
                    aComHSfM.AddArgs("InMulTieP="+mPhProj.DPMulTieP().DirOut());
                    continue;
                }

            }
            if (starts_with(mArgv[aKP], "OutOri")) // already added
                continue;
            if (starts_with(mArgv[aKP], "VirTP")) // already added
                continue;
            if (starts_with(mArgv[aKP], "Exe"))
                continue;


            aComHSfM.AddArgs(mArgv[aKP]);
        }
        aLPCS.push_back(aComHSfM);
    }


    if (mExe)
        ExeComSerial(aLPCS, true);
    else
    {
        for (auto aCmd : aLPCS)
            StdOut() << aCmd.Com() << std::endl;
    }

    return EXIT_SUCCESS;
}

/* ====================================================== */
/*               OriPoseEstimRelGlob                       */
/* ====================================================== */

tMMVII_UnikPApli Alloc_OriRelPipeline(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
    return tMMVII_UnikPApli(new cAppli_OriRelPipeline(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_OriPoseEstimRel
    (
        "OriPoseEstimRelGlob",
        Alloc_OriRelPipeline,
        "Global orientation using relative motions (pairs, triplets)",
        {eApF::Ori},
        {eApDT::TieP},
        {eApDT::Orient},
        __FILE__
    );

}


