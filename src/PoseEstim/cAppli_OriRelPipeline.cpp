#include "MMVII_DeclareAllCmd.h"
#include "MMVII_Sensor.h"
#include "MMVII_PoseRel.h"

namespace MMVII
{

const std::vector<std::string> & Vec_NP_VirTP()
{
    static    std::vector<std::string> aRes {NP_VirTP,NP_DistribVirTP} ;

    return aRes;
}

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
        eDistrVirTPs              mDistribVirTPs;
        bool                      mExe;
        bool                      mSelTriCov;
        bool                      mCalledByEstimFocal;
        eStepOriRelPipeline       mStepEnd;
};

cAppli_OriRelPipeline::cAppli_OriRelPipeline(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
    cMMVII_Appli  (aVArgs,aSpec),
    mPhProj (*this),
    mDistribVirTPs(eDistrVirTPs::e5Pts),
    mExe(true),
    mSelTriCov(false),
    mCalledByEstimFocal (false),
    mStepEnd            (eStepOriRelPipeline::eOriAll)
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
           <<  AOpt2007(mExe,"Exe","Execute pipeline",{eTA2007::HDV})
           <<  AOpt2007(mSelTriCov,"SelTriCov","Select triplets by tie-points coverage and B/H (OriPoseSelecAllTriplets)",{eTA2007::HDV})
           << cHeaderSectionArg("Input for tie points (at least one must be provided)")
           << mPhProj.AddAllDpTieP(anArgOpt)

           << cHeaderSectionArg("Export virtual points")
           <<  mPhProj.DPMulTieP().ArgDirOutOpt(NP_VirTP,"Output folder for virtual tie points")
           << AOpt2007(mDistribVirTPs,NP_DistribVirTP,"Distribution of virtual tie points",{eTA2007::HDV})

           << cHeaderSectionArg("Run global orientation")
           <<  mPhProj.DPOrient().ArgDirOutOpt("","Output Ori: if set, the global orientation is computed")

            << cHeaderSectionArg("For Call by Estim Focal")
            << AOpt2007(mCalledByEstimFocal,NP_CalledByEstimFocal,"Signal a call by estim focal",{eTA2007::HDV,eTA2007::ExternalCall})
            << AOpt2007(mStepEnd,NP_StepEndEstimFocal,"Step where we stop",{eTA2007::HDV,eTA2007::ExternalCall})
        ;
}

int cAppli_OriRelPipeline::Exe()
{
    mPhProj.FinishInit();

    bool  PreviousBehaviour = false;

    mPhProj.AssertAtLeastOneTiePModeIsInit();
    /*
    MMVII_INTERNAL_ASSERT_always((mPhProj.DPTieP().DirInIsInit() ||
                                  mPhProj.DPMulTieP().DirInIsInit() ||
                                  mPhProj.DPGndPt2D().DirInIsInit()),"Input tie-points not provided");*/



    std::list<cParamCallSys>  aLPCS;

    // ===============================================
    // select pairs

    if (PreviousBehaviour)
    {

        std::string aTieP = mPhProj.DPTieP().DirInIsInit() ? ("InTieP="+mPhProj.DPTieP().DirIn()) :
                        (mPhProj.DPMulTieP().DirInIsInit() ? ("InMulTieP="+mPhProj.DPMulTieP().DirIn()) :
                                                             ("InObjMesInstr="+mPhProj.DPGndPt2D().DirIn()) );
        cParamCallSys aComSelPairs(cMMVII_Appli::FullBin(),
                                   "OriPoseSelecAllPAir",
                                   mPatIm,
                                   mPhProj.DPOriRel().DirOut(),
                                   aTieP);
        aLPCS.push_back(aComSelPairs);
    }
    else
    {
        cParamCallSys aComSelPairs(TheSpec_SelectPairOriRel,true);

        aComSelPairs.AddMandArgsOfAppli(*this,{0,1}); // Add 2 first args
        aComSelPairs.AddMMVIIAllTiepInInit(mPhProj);   // all the parameters relative to TieP
        if (mCalledByEstimFocal)
           aComSelPairs.AddNoMsg() ;
        aLPCS.push_back(aComSelPairs);
    }


    // ===============================================
    // compute relative orientations between pairs
    if (PreviousBehaviour)
    {

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
    }
    else
    {
        cParamCallSys aComEstimRelPair(TheSpec_OriRelAllPairs,true);

        // aComEstimRelPair.AddAllMandArgsOfAppli(*this);
        aComEstimRelPair.AddMandArgsOfAppli(*this,{1,2});
        aComEstimRelPair.AddMMVIIAllTiepInInit(mPhProj);   // all the parameters relative to TieP
        aComEstimRelPair.AddArgsOfAppli(*this,Vec_NP_VirTP());  // All the parmeter relative to Vir Tiep
        if (mCalledByEstimFocal)
           aComEstimRelPair.AddNoMsg() ;

        if (int(mStepEnd)>= (int) eStepOriRelPipeline::eOriPairs)
           aLPCS.push_back(aComEstimRelPair);
    }

    // ===============================================
    // select triplets by coverage and B/H
    if (mSelTriCov && (int(mStepEnd) >= (int) eStepOriRelPipeline::eOriTriplets))
    {
        // the pair selection below rewrites the cycle triplets : the backup of a previous run is stale
        if (mExe)
        {
            std::string aPostCyc = "Cyc." + GlobTaggedNameDefSerial();
            for (const auto & aNameIm : VectMainSet(0))
            {
                std::string aNameCyc = mPhProj.OriRel_NameAllTripletsOf1Image(aNameIm,false,aPostCyc);
                if (ExistFile(aNameCyc))
                    RemoveFile(aNameCyc,SVP::No);
            }
        }

        cParamCallSys aComSelTri(TheSpec_SelecAllTriplets,true);

        aComSelTri.AddMandArgsOfAppli(*this,{1,2});
        aComSelTri.AddMMVIIAllTiepInInit(mPhProj);   // all the parameters relative to TieP
        if (mCalledByEstimFocal)
            aComSelTri.AddNoMsg();
        aLPCS.push_back(aComSelTri);
    }

    // ===============================================
    // compute relative orientations within triplets
    if (PreviousBehaviour)
    {
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
    }
    else
    {
        cParamCallSys aComEstimRelTri(TheSpec_OriRelAllTriplets,true);

        // aComEstimRelPair.AddAllMandArgsOfAppli(*this);
        aComEstimRelTri.AddMandArgsOfAppli(*this,{1,2});
        aComEstimRelTri.AddMMVIIAllTiepInInit(mPhProj);   // all the parameters relative to TieP
        aComEstimRelTri.AddArgsOfAppli(*this,Vec_NP_VirTP());
        if (mCalledByEstimFocal)
           aComEstimRelTri.AddNoMsg() ;

        if (int(mStepEnd)>= (int) eStepOriRelPipeline::eOriTriplets)
            aLPCS.push_back(aComEstimRelTri);
    }


    // ===============================================
    // global hierarchical SFM
    if (mPhProj.DPOrient().DirOutIsInit())
    {
        if (PreviousBehaviour)
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
                    // MPD=>ER : la il y a un pb, pk si avait utilise des tiep "normaux" ou des cible, ca passe pas
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
                if (starts_with(mArgv[aKP], "SelTriCov"))
                    continue;


                aComHSfM.AddArgs(mArgv[aKP]);
            }
            aLPCS.push_back(aComHSfM);
        }
        else
        {
            cParamCallSys aComHSfM(TheSpec_HierarchSfm,true);

            // the vir tie points are  given as  mPhProj.DPMulTieP().ArgDirOutOpt()
            bool UseVirtTieP = mPhProj.DPMulTieP().DirOutIsInit();

            aComHSfM.AddMandArgsOfAppli(*this,{0,1,2});     // Add Pattern +
            aComHSfM.AddArgs(mPhProj.DPOrient().DirOut());  // The optionnal dirout is a mandatoy param of SFM

            if (UseVirtTieP)
            {
                // Extract name  "InMulTieP" using standar interface
                std::string aNameParam = cDirsPhProj::DefNameArgOpt(eTA2007::MulTieP,true);
                aComHSfM.AddMMVIIArgsOpt(aNameParam,mPhProj.DPMulTieP().DirOut());
            }
             else
            {
                // as usual : add all the parameters relative to TieP
               aComHSfM.AddMMVIIAllTiepInInit(mPhProj);
            }

            if (mCalledByEstimFocal)
               aComHSfM.AddNoMsg() ;

            if (int(mStepEnd)>= (int) eStepOriRelPipeline::eOriAll)
               aLPCS.push_back(aComHSfM);
        }
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


