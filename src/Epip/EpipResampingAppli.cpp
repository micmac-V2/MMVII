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
    const cInterpolator1D* anInterp = cDiffInterpolator1D::AllocFromNames(mInterpol);


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

    const auto* Im1 = ReadIm2DGen(mNameIm1);
    const auto* Im2 = ReadIm2DGen(mNameIm2);

    StdOut() << "Interpolator: " << anInterp->VNames() << ", Kernel Size: " << anInterp->SzKernel() << std::endl;

    auto aResult = generateEpipolarImages(aEpipModel, Im1, Im2, *anInterp, mMargin);


    // TODOCM: Make name generation accessible for other apps
    // TODOCM: Make sure image extension is present (and not doubled ?) ! (i.e. .tif)
    auto anEpip1Name = aOutDir + replaceFirstOccurrence(replaceFirstOccurrence(mOutNamePat,"%1",mNameIm1),"%2",mNameIm2);
    auto anEpip2Name = aOutDir + replaceFirstOccurrence(replaceFirstOccurrence(mOutNamePat,"%1",mNameIm2),"%2",mNameIm1);
    aResult.im1_rect->ToFile(anEpip1Name);
    aResult.im2_rect->ToFile(anEpip2Name);
    StdOut() << "Image1: " << anEpip1Name << std::endl;
    StdOut() << "Image2: " << anEpip2Name << std::endl;


    delete aResult.im1_rect;
    delete aResult.im2_rect;
    delete anInterp;
    delete Im1;
    delete Im2;
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

