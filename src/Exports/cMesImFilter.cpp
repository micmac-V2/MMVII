#include "cMesImFilter.h"
//#include "MMVII_Sensor.h"
#include "MMVII_ImageMorphoMath.h"
#include "MMVII_Interpolators.h"
#include "MMVII_PCSens.h"
#include "cMMVII_Appli.h"
#include "MMVII_Geom3D.h"
#include "MMVII_Matrix.h"
#include "MMVII_Tpl_ElemStrToVal.h"


namespace MMVII
{

struct cImFiltParam
{
    std::string ImFilt;
    bool Reverse;
    ARG2007_STRUCT_FIELDS(
        ImFilt, FieldSem({eTA2007::AddCom,"Only affect image with name that match pattern"}),
        Reverse, FieldSem({eTA2007::AddCom,"Reverse pattern"})
    )
};

class cAppli_MesImFilter : public cMMVII_Appli//heritage de cMMVII_Appli
{
public:
    //public method/attributes declaration
    cAppli_MesImFilter(const std::vector<std::string>& aVArgs,
                     const cSpecMMVII_Appli& aSpec);
private:
    //private method/attributes declaration
    int Exe() override;
    cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
    cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;
    cPhotogrammetricProject mPhProj;
    std::string mSpecImIn;
    std::string mMesFilt;

    bool mShow;
    std::string mMesOut;
    cImFiltParam mImFilt;
};

cCollecSpecArg2007& cAppli_MesImFilter::ArgObl(cCollecSpecArg2007& anArgObl)
{
    return anArgObl
           << Arg2007(mSpecImIn, "Pattern/file of images", {{eTA2007::MPatFile,"0"}, {eTA2007::FileDirProj}})
           << mPhProj.DPGndPt2D().ArgDirInMand("Input image measurements")
        ;
}

cCollecSpecArg2007 & cAppli_MesImFilter::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    return anArgOpt
           << mPhProj.DPGndPt2D().ArgDirOutOpt(mMesOut, "Ouptut image measurements", true)
           << AOpt2007(mImFilt, "ImFilt", "Filter only specific images", {eTA2007::HDV, eTA2007::AddCom})
           << AOpt2007(mShow, "Show", "Show execution details", {eTA2007::HDV})
        ;
}

cAppli_MesImFilter::cAppli_MesImFilter(const std::vector<std::string>& aVArgs,
                                   const cSpecMMVII_Appli& aSpec):
    cMMVII_Appli(aVArgs, aSpec),
    mPhProj (*this),
    mImFilt ({"*", false})
{
    //constructor does nothing
}

int cAppli_MesImFilter::Exe()
{
    mPhProj.FinishInit();
    std::vector<std::string> aVIm = VectMainSet(0);
    return EXIT_SUCCESS;
}

tMMVII_UnikPApli Alloc_MesImFilter(const std::vector<std::string> & aVArgs,
                                 const cSpecMMVII_Appli & aSpec)
{
    return tMMVII_UnikPApli(new cAppli_MesImFilter(aVArgs, aSpec));
}

cSpecMMVII_Appli TheSpec_MesImFilter
    (
        "MesImFilter",
        Alloc_MesImFilter,
        "generate target 3D description from poses & images measurements",
        //metadonnees
        {eApF::Ori,eApF::GCP},//features
        {eApDT::ObjCoordWorld, eApDT::ObjMesInstr},//inputs
        {eApDT::Console},//output
        __FILE__
        );
}

