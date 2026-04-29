#include "cMMVII_Appli.h"
#include "MMVII_Sensor.h"
#include "cEpipolarRectification.h"
#include "MMVII_Interpolators.h"
#include <vector>

/**
   \file EpipGeom.cpp


 */


namespace MMVII
{


class cAppli_EpipResampling : public cMMVII_Appli
{
public :

    cAppli_EpipResampling(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli &);
    int Exe() override;
    cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
    cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;

private :
    cPhotogrammetricProject  mPhProj;
    std::string  mNameIm1;
    std::string  mNameIm2;
    int mDegree = 5;
    int mDegreeInv = mDegree + 4;
    int mNbByXY = 100;
    int mNbByZ = 3;
    double mEpsMarginRel = 0.05;
    int mMargin = 2;
    std::string mOutDir;
    std::string mOutNamePat = "Epip_%1_%2.tif";
    std::vector<std::string> mInterpol = {"Cubic","-0.5"};
};

cAppli_EpipResampling::cAppli_EpipResampling (
    const std::vector<std::string> &  aVArgs,
    const cSpecMMVII_Appli & aSpec
    )
    : cMMVII_Appli  (aVArgs,aSpec)
    , mPhProj       (*this)
{
}


// TODOCM: Separate Geom calulation & resampling
// TODOCM: Serialisation classes EpipMap EpipModel
// TODOCM: Creer les RPC des images Epipolaires
// TODOCM: X Steps /= Y Steps. Steps or pixels  => degre liberté * 10 , bonne répartition, Nb min en X et Y.

// TODOCM: Gestion grosses images : daller ... Cache pour bout d'images ?

// TODOCM: adapter/utiliser le resample de la cDataIm2D




int cAppli_EpipResampling::Exe()
{
    mPhProj.FinishInit();

    if (! IsInit(&mDegreeInv))
        mDegreeInv = mDegree + 4;

    std::string aOutDir = mPhProj.DirVisuAppli();
    if (IsInit(&mOutDir))
    {
        aOutDir = mOutDir;
    }
    if (! aOutDir.empty())
    {
        aOutDir += "/";
    }
    CreateDirectories(aOutDir,false);
    const cInterpolator1D* aInterp = cDiffInterpolator1D::AllocFromNames(mInterpol);


    // TODOCM: Enlever margin ? Mieux le définir ?
    // TODOCM: Test d'epipolarisabilite ...


    const cSensorImage *  aSI1 =  mPhProj.ReadSensor(FileOfPath(mNameIm1,false /* Ok Not Exist*/),true/*DelAuto*/,false /* Not SVP*/);
    if (! aSI1->HasIntervalZ())
    {
        MMVII_UserError(eTyUEr::eOpenFile,"Image 1 is not from a RPC sensor");
    }
    const cSensorImage *  aSI2 =  mPhProj.ReadSensor(FileOfPath(mNameIm2,false /* Ok Not Exist*/),true/*DelAuto*/,false /* Not SVP*/);
    if (! aSI2->HasIntervalZ())
    {
        MMVII_UserError(eTyUEr::eOpenFile,"Imag,e 2 is not from a RPC sensor");
    }

    auto aDIm1 = cDataFileIm2D::Create(mNameIm1,eForceGray::No);
    auto aDIm2 = cDataFileIm2D::Create(mNameIm2,eForceGray::No);
    StdOut() <<  "Image1=" <<  mNameIm1;
    StdOut() << " " << aDIm1.Sz() << " " << ToStr(aDIm1.Type()) << " " << aDIm1.NbChannel() << " chan" << std::endl;
    StdOut() <<  "Image2=" <<  mNameIm2;
    StdOut() << " " << aDIm2.Sz() << " "  << ToStr(aDIm2.Type()) << " " << aDIm2.NbChannel() << " chan" << std::endl;

    auto aParams = cEpipolarRectification::cParams{mDegree,mDegreeInv,mNbByXY,mNbByZ,mEpsMarginRel};
    auto aRectifier = cEpipolarRectification(*aSI1, *aSI2, aParams);
    auto aEpipModel = aRectifier.Compute();

    const auto* aIm1 = ReadIm2DGen(mNameIm1);
    const auto* aIm2 = ReadIm2DGen(mNameIm2);

    aEpipModel.ComputeCommonFraming(aIm1->Sz(),aIm2->Sz());
    std::cout << "[EpipolarFraming]\n"
              << "  Image 1: rectified X=[" << aEpipModel.Frame1.xMin_rect << ", " << aEpipModel.Frame1.xMax_rect
              << "]  Y=[" << aEpipModel.Frame1.yMin_rect << ", " << aEpipModel.Frame1.yMax_rect << "]\n"
              << "  Image 2: rectified X=[" << aEpipModel.Frame2.xMin_rect << ", " << aEpipModel.Frame2.xMax_rect
              << "]  Y=[" << aEpipModel.Frame2.yMin_rect << ", " << aEpipModel.Frame2.yMax_rect << "]\n"
              << "  Common Y range : [" << aEpipModel.Frame1.yMin_common << ", " << aEpipModel.Frame1.yMax_common << "]\n"
              << "  Common outSy   : " << aEpipModel.Frame1.outSy << " px\n"
              << "  Image 1 outSx  : " << aEpipModel.Frame1.outSx << " px\n"
              << "  Image 2 outSx  : " << aEpipModel.Frame2.outSx << " px\n";


    StdOut() << "Interpolator: " << aInterp->VNames() << ", Kernel Size: " << aInterp->SzKernel() << std::endl;

    std::cout << "[EpipolarResample] Resampling image 1 ("
              << aEpipModel.Frame1.outSx << "x" << aEpipModel.Frame1.outSy << ")...\n";
    auto aIm1Rectif = aEpipModel.Resample1(aIm1,*aInterp);

    std::cout << "[EpipolarResample] Resampling image 2 ("
              << aEpipModel.Frame2.outSx << "x" << aEpipModel.Frame2.outSy << ")...\n";
    auto aIm2Rectif = aEpipModel.Resample2(aIm2,*aInterp);


    // TODOCM: Make name generation accessible for other apps
    // TODOCM: Make sure image extension is present (and not doubled ?) ! (i.e. .tif)
    auto anEpip1Name = aOutDir + replaceFirstOccurrence(replaceFirstOccurrence(mOutNamePat,"%1",mNameIm1),"%2",mNameIm2);
    auto anEpip2Name = aOutDir + replaceFirstOccurrence(replaceFirstOccurrence(mOutNamePat,"%1",mNameIm2),"%2",mNameIm1);
    aIm1Rectif->ToFile(anEpip1Name);
    aIm2Rectif->ToFile(anEpip2Name);
    StdOut() << "Image1: " << anEpip1Name << std::endl;
    StdOut() << "Image2: " << anEpip2Name << std::endl;

    auto a = aIm1->AllocReSampleGen(*aInterp,aEpipModel.EpipMap1());
    auto b = aIm2->AllocReSampleGen(*aInterp,aEpipModel.EpipMap2());

    a.second->ToFile("Joe1");
    b.second->ToFile("Joe2");

    delete aIm1Rectif;
    delete aIm2Rectif;
    delete aInterp;
    delete aIm1;
    delete aIm2;
    return EXIT_SUCCESS;
}


cCollecSpecArg2007 & cAppli_EpipResampling::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return anArgObl
          << Arg2007(mNameIm1,"name first image",{eTA2007::FileImage})
          << Arg2007(mNameIm2,"name second image",{eTA2007::FileImage})
          << mPhProj.DPOrient().ArgDirInMand()
        ;
}

extern cSpecMMVII_Appli  TheSpec_EpipResampling;  // Forward declaration

cCollecSpecArg2007 & cAppli_EpipResampling::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{

    return anArgOpt
           << AOpt2007(mDegree,"Degree","Poly degree",{eTA2007::HDV})
           << AOpt2007(mDegreeInv,"DegreeInv","Inv Poly degree",{eTA2007::HDV})
           << AOpt2007(mNbByXY,"XYSteps","Nb XY steps",{eTA2007::HDV})
           << AOpt2007(mNbByZ,"ZSteps","Nb Z steps",{eTA2007::HDV})
           << AOpt2007(mEpsMarginRel,"MarginRel","Relative margin for H-compatible points grid (X,Y,Z)",{eTA2007::HDV})
           << AOpt2007(mMargin,"Margin","Output image margin (black contour)",{eTA2007::HDV})
           << AOpt2007(mOutDir,"OutDir","Output directory (Default: VISU/" + TheSpec_EpipResampling.Name()+")")
           << AOpt2007(mOutNamePat,"OutNamer","Output name pattern", {eTA2007::HDV})
           << AOpt2007(mInterpol,"Interpol","Interpolator", {eTA2007::HDV})
        ;
}



/* ==================================================== */

tMMVII_UnikPApli Alloc_EpipResampling(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
    return tMMVII_UnikPApli(new cAppli_EpipResampling(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_EpipResampling
    (
        "EpipResampling",
        Alloc_EpipResampling,
        "Epipolar geometry of two images",
        {eApF::Ori},
        {eApDT::Orient},
        {eApDT::Orient},
        __FILE__
        );


}; // MMVII

