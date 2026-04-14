#include "cCodedTargetDescribe.h"

namespace MMVII
{

    cCollecSpecArg2007& cAppli_CodedTargetDescribe::ArgObl(cCollecSpecArg2007& anArgObl)
    {
        return anArgObl
               << Arg2007(mSpecImIn, "Pattern/file of images", {{eTA2007::MPatFile,"0"}, {eTA2007::FileDirProj}})
               << mPhProj.DPOrient().ArgDirInMand()
            ;
    }

    cCollecSpecArg2007 & cAppli_CodedTargetDescribe::ArgOpt(cCollecSpecArg2007 & anArgOpt)
    {
        return anArgOpt
               << AOpt2007(mShow,"Show","show some useful details", {eTA2007::HDV})//hdv = has default value
            ;
    }

    cAppli_CodedTargetDescribe::cAppli_CodedTargetDescribe(const std::vector<std::string>& aVArgs,
                                                           const cSpecMMVII_Appli& aSpec):
        cMMVII_Appli(aVArgs, aSpec),
        mPhProj(*this)
    {
        //constructor does nothing
    }

    int cAppli_CodedTargetDescribe::Exe()
    {
        //----- [0] Load project primitives
        mPhProj.FinishInit();
        std::vector<std::string> aVIm = VectMainSet(0);
        //----- [1] Load targets measures
        mPhProj.LoadGCP3D();
        for (const std::string& aIm : aVIm)
        {
            //cSensorCamPC* aCam = mPhProj.ReadCamPC(aIm, true);
            cSetMesPtOf1Im aSetImMes = mPhProj.LoadMeasureIm(aIm);
            std::vector<cSaveExtrEllipe> aVEllipsesExtrinsics;
            ReadFromFile(aVEllipsesExtrinsics, cSaveExtrEllipe::NameFile(mPhProj, aSetImMes, true));
        }
        StdOut() << "OK";
        return EXIT_SUCCESS;
    }

    tMMVII_UnikPApli Alloc_CodedTargetDescribe(const std::vector<std::string> & aVArgs,
                                                  const cSpecMMVII_Appli & aSpec)
    {
        return tMMVII_UnikPApli(new cAppli_CodedTargetDescribe(aVArgs, aSpec));
    }

    cSpecMMVII_Appli TheSpec_CodedTargetDescribe
        (
            "CodedTargetDescribe",
            Alloc_CodedTargetDescribe,
            "generate target 3D description from poses & images measurements",
            //metadonnees
            {eApF::Ori,eApF::GCP},//features
            {eApDT::ObjCoordWorld, eApDT::ObjMesInstr},//inputs
            {eApDT::Console},//output
            __FILE__
        );


    }
