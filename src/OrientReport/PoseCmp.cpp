#include "MMVII_Ptxd.h"
#include "cMMVII_Appli.h"
#include "MMVII_Geom3D.h"
#include "MMVII_PCSens.h"
#include "MMVII_Tpl_Images.h"
#include "MMVII_StaticLidar.h"

/**
   \file GCPQuality.cpp


 */

namespace MMVII
{

class cAppli_PoseCmp : public cMMVII_Appli
{
public :

    cAppli_PoseCmp(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli &);
    int Exe() override;
    cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
    cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;

private :
    bool AcceptEmptySet(int aK) const override {return ((aK==0)&&(mSpecImIn==MMVII_NONE));}
    cPhotogrammetricProject  mPhProj;
    std::vector<int>         mPropStat;
    std::string              mSpecImIn;   ///  Pattern of xml file
    std::vector<std::string> mSetNames;

    std::string              mOri2;
    std::vector<std::string> mPatBand;
    bool                     mDoRel;
    bool                     mVerbose; /// show each difference

    // std::string              mPrefixCSV;
    // std::string              mPrefixCSVIma;

    void MakeStatByImage();
};

cAppli_PoseCmp::cAppli_PoseCmp
    (
        const std::vector<std::string> &  aVArgs,
        const cSpecMMVII_Appli & aSpec
        ) :
    cMMVII_Appli  (aVArgs,aSpec),
    mPhProj       (*this),
    mPropStat     ({50,75}),
    mDoRel        (false)
{
}

cCollecSpecArg2007 & cAppli_PoseCmp::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return     anArgObl
           << Arg2007(mSpecImIn,"Pattern/file for images",{{eTA2007::MPatFile,"0"},{eTA2007::FileDirProj}})
           << mPhProj.DPOrient().ArgDirInMand()
           << Arg2007(mOri2,"Second orientation folder",{{eTA2007::Input},{eTA2007::Orient}})
        ;
}


cCollecSpecArg2007 & cAppli_PoseCmp::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{

    return   anArgOpt
           << AOpt2007(mPropStat,"Perc","Percentil for stat exp",{eTA2007::HDV})
           << AOpt2007(mPatBand,"PatBand","Pattern for band [Patter,Subts]")
           << AOpt2007(mDoRel,"DoRel","Do relative computaion between consecutive images")
           << AOpt2007(mVerbose, "Verbose"  ,"Do show each difference",{eTA2007::HDV})
        ;
}



int cAppli_PoseCmp::Exe()
{
    mPhProj.FinishInit();

    cDirsPhProj &  aDirOri2 = *mPhProj.NewDPIn(eTA2007::Orient,mOri2);

    MMVII_INTERNAL_ASSERT_tiny(IsDirectory(mPhProj.DPOrient().FullDirIn()),
                               "Ori " + mPhProj.DPOrient().FullDirIn() + " does not exist");

    MMVII_INTERNAL_ASSERT_tiny(IsDirectory(aDirOri2.FullDirIn()),
                               "Ori " + mOri2 + " does not exist");

    std::vector<cSensorCamPC*> aVCam1;
    std::vector<cSensorCamPC*> aVCam2;

    std::string aLastBand = "Xy#@Z-4lj";

    cPt3dr aLastWPK=cPt3dr::Dummy();
    cWeightAv<tREAL8>  aAvgDif_Ori;
    cWeightAv<tREAL8>  aAvgDif_Center;
    cWeightAv<tREAL8>  aAvgRelDif_Ori;
    cWeightAv<tREAL8>  aAvgBandRelDif;

    for (const auto & aNameIm : VectMainSet(0))
    {
        cSensorCamPC * aCam1 = mPhProj.ReadCamPC(aNameIm,true,SVP::Yes);
        cSensorCamPC * aCam2 = mPhProj.ReadCamPC(aDirOri2,aNameIm,true,SVP::Yes);
        if (aCam1==nullptr)
            aCam1 =  mPhProj.ReadStaticLidar(aNameIm,true,SVP::Yes,false);
        if (aCam2==nullptr)
            aCam2 =  mPhProj.ReadStaticLidar(aDirOri2,aNameIm,true,SVP::Yes,false);
        if (!aCam1 || !aCam2)
            continue;

        aVCam1.push_back(aCam1);
        aVCam2.push_back(aCam2);

        auto aP2In1 = aCam1->RelativePose(*aCam2);
        cPt3dr aWPK = aP2In1.Rot().ToWPK() ;
        aAvgDif_Ori.Add(1.0,Norm2(aWPK));
        aAvgDif_Center.Add(1.0,Norm2(aP2In1.Tr()));
        if (mVerbose)
            StdOut()<<aNameIm<<": Ori="<<Norm2(aWPK)<<", Center="<<Norm2(aP2In1.Tr())<<"\n";

        if (aLastWPK.IsValid())
            aAvgRelDif_Ori.Add(1.0,Norm2(aWPK-aLastWPK));

        if (IsInit(&mPatBand))
        {
            std::string aNameBand = ReplacePattern(mPatBand.at(0),mPatBand.at(1),aNameIm);
            bool aNewBand = (aNameBand!=aLastBand);
            aLastBand = aNameBand;

            if (aNewBand)
                StdOut() << "NB=" << aNameIm << std::endl;
            else
                aAvgBandRelDif.Add(1.0,Norm2(aWPK-aLastWPK));
        }

        // StdOut() <<  "Im=" << aName << " Tr="  << aP2In1.Tr()  << " WPK=" << aP2In1.Rot().ToWPK() << std::endl;
        aLastWPK = aWPK;
    }

    // mPrefixCSV =  "_Ori-"+  mPhProj.DPOrient().DirIn() +  "_Mes-"+  mPhProj.DPMulTieP().DirIn() ;
    //
    StdOut() << mPhProj.DPOrient().DirIn() << " " << mOri2 << ": ";
    StdOut() << "AVG DIFF, Ori=" << aAvgDif_Ori.Average() << " Center=" << aAvgDif_Center.Average() << std::endl;
    if (mDoRel)
        StdOut() << "AVG REL DIFF=" << aAvgRelDif_Ori.Average() << std::endl;
    if (aAvgBandRelDif.SW() > 0)
        StdOut() << "AVG BAND REL DIFF=" << aAvgBandRelDif.Average() << std::endl;

    return EXIT_SUCCESS;
}


/* ==================================================== */

tMMVII_UnikPApli Alloc_PoseCmp(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
    return tMMVII_UnikPApli(new cAppli_PoseCmp(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_PoseCmpReport
    (
        "ReportPoseCmp",
        Alloc_PoseCmp,
        "Reports on pose comparison",
        {eApF::Ori},
        {eApDT::Orient},
        {eApDT::Csv},
        __FILE__
    );


}; // MMVII

