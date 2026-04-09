//#include "MMVII_Sensor.h"
#include "CodedTarget.h"
#include "MMVII_ImageMorphoMath.h"
#include "MMVII_Interpolators.h"
#include "MMVII_PCSens.h"
#include "cMMVII_Appli.h"
#include "MMVII_Geom3D.h"
#include "MMVII_Matrix.h"

namespace MMVII
{

    class cAppli_CodedTargetDescribe : public cMMVII_Appli//heritage de cMMVII_Appli
    {
        public:
            //public method/attributes declaration
            cAppli_CodedTargetDescribe(const std::vector<std::string>& aVArgs,
                                       const cSpecMMVII_Appli& aSpec);
        private:
            //MMVII mandatory/usual stuff
            int Exe() override;
            cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
            cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;
            cPhotogrammetricProject mPhProj;
            std::string mSpecImIn;
            bool mShow;
            //------
            void establishPrerequisites();
    };

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
