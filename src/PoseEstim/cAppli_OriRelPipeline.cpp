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
        bool                      mExe;
};

cAppli_OriRelPipeline::cAppli_OriRelPipeline(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
    cMMVII_Appli  (aVArgs,aSpec),
    mPhProj (*this),
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
           <<  mPhProj.DPTieP().ArgDirInOpt()
           <<  mPhProj.DPGndPt2D().ArgDirInOpt()
           <<  mPhProj.DPMulTieP().ArgDirInOpt()
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




    // compute relative orientations between pairs
    cParamCallSys aComEstimRelPair(cMMVII_Appli::FullBin(),
                                   "OriPoseEstimRelAllPairs");

    for (size_t aKP=3; aKP<mArgv.size() ; aKP++)
    {
        if (!  starts_with(mArgv[aKP], "VirTP") && ! starts_with(mArgv[aKP], "Exe")) //omit virtual points for pairs
            aComEstimRelPair.AddArgs(mArgv[aKP]);
    }
    aLPCS.push_back(aComEstimRelPair);





    // compute relative orientations within triplets
    cParamCallSys aComEstimRelTri(cMMVII_Appli::FullBin(),
                                   "OriPoseEstimRelAllTriplets");

    for (size_t aKP=3; aKP<mArgv.size() ; aKP++)
    {
        if (! starts_with(mArgv[aKP], "Exe"))
            aComEstimRelTri.AddArgs(mArgv[aKP]);
    }
    aLPCS.push_back(aComEstimRelTri);




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
/*               OriPoseEstimRel                       */
/* ====================================================== */

tMMVII_UnikPApli Alloc_OriRelPipeline(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
    return tMMVII_UnikPApli(new cAppli_OriRelPipeline(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_OriPoseEstimRel
    (
        "OriPoseEstimRel",
        Alloc_OriRelPipeline,
        "Estimate relative orientations of all images (pairs and triplets)",
        {eApF::Ori},
        {eApDT::TieP},
        {eApDT::Orient},
        __FILE__
    );

}


