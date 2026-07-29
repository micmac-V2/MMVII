#include "cMMVII_Appli.h"
#include "MMVII_Sensor.h"
#include "cEpipolarRectification.h"
#include "MMVII_Interpolators.h"
#include "MMVII_CodeTiming.h"
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
    int mMargin = 2;
    std::string mOutDir;
    std::string mOutNamePat = "Epip_%1_%2.tif";
    std::vector<std::string> mInterpol = {"Cubic","-0.5"};
    eEpipFrm mFrame = eEpipFrm::eIntersect;
//    cEpipolarModel::eFramingType mFrame = cEpipolarModel::eFramingType::INTERSECT;
};

cAppli_EpipResampling::cAppli_EpipResampling (
    const std::vector<std::string> &  aVArgs,
    const cSpecMMVII_Appli & aSpec
    )
    : cMMVII_Appli  (aVArgs,aSpec)
    , mPhProj       (*this)
{
}


// TODOCM: Serialisation classes EpipMap EpipModel ?
// TODOCM: X Steps /= Y Steps. Steps or pixels  => degre liberté * 10 , bonne répartition, Nb min en X et Y.
// TODOCM: Gestion grosses images : daller ... Cache pour bout d'images ?
// TODOCM: Enlever margin ? Mieux le définir ?
// TODOCM: Test d'epipolarisabilite ...
// TODOCM: Dans GenerateData, ajouter un point avec Z aléatoire
// TODOCM: Dans GenerateData creer plus de point et utiliser un point sur 2 pour le calcul des poly et les autres pour la verif
// TODOCM: Prendre en compte le fait que les 2 images n'ont pas forcement le même intervalle Z
// TODOCM: Faire des benchs en utilisant recalcul des RPCs
// TODOCM: Make output image name generation accessible for other apps
// TODOCM: Make sure image  file extension is present ? (.tif)
// TODOCM: Changer nom fichier RPC (RPC-xxx.xml ?)
// TODOCM: Vraiment tester tabulation pour EpipModel


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
    CreateDirectories(aOutDir);
    const cInterpolator1D* aInterp = cDiffInterpolator1D::AllocFromNames(mInterpol);


    const cSensorImage *  aSI1 =  mPhProj.ReadSensor(FileOfPath(mNameIm1,false /* Ok Not Exist*/),true/*DelAuto*/,false /* Not SVP*/);
    if (! aSI1->HasIntervalZ())
    {
        MMVII_UserError(eTyUEr::eOpenFile,"Image 1 has no Z validity interval");
    }

    const cSensorImage *  aSI2 =  mPhProj.ReadSensor(FileOfPath(mNameIm2,false /* Ok Not Exist*/),true/*DelAuto*/,false /* Not SVP*/);
    if (! aSI2->HasIntervalZ())
    {
        MMVII_UserError(eTyUEr::eOpenFile,"Image 2 has no Z validity interval");
    }

    // Create early to have possible error reported before doing long computations
    auto aDIm1 = cDataFileIm2D::Create(mNameIm1,eForceGray::No);
    auto aDIm2 = cDataFileIm2D::Create(mNameIm2,eForceGray::No);

    StdOut() <<  "Image_1: " <<  mNameIm1;
    StdOut() << " " << aDIm1.Sz() << " " << ToStr(aDIm1.Type()) << " " << aDIm1.NbChannel() << " chan" << std::endl;
    StdOut() <<  "Image_2: " <<  mNameIm2;
    StdOut() << " " << aDIm2.Sz() << " "  << ToStr(aDIm2.Type()) << " " << aDIm2.NbChannel() << " chan" << std::endl;

    StdOut() << "Degree: " << mDegree << ", DegreeInv: " << mDegreeInv << std::endl;
    StdOut() << "NbByXY: " << mNbByXY << ", NbByZ: " << mNbByZ << std::endl;
    StdOut() << "Frame: " << ToStr(mFrame) << ", Margin: " << mMargin << std::endl;

    auto aParams = cEpipolarRectification::cParams{mDegree,mDegreeInv,mNbByXY,mNbByZ,mFrame,mMargin};
    auto aRectifier = cEpipolarRectification(*aSI1, *aSI2, aParams);
    auto aEpipModel = aRectifier.Compute();

    StdOut() << "Nb Pairs 1->2 : " << aRectifier.NbPairs12() << std::endl;
    StdOut() << "Nb Pairs 2->1 : " << aRectifier.NbPairs21() << std::endl;
    StdOut() << "V1,V2 var : " << aRectifier.V1V2Var() << std::endl;
    StdOut() << "W1 var : " << aRectifier.W1Var() << std::endl;
    StdOut() << "W2 var : " << aRectifier.W2Var() << std::endl;


    const auto& anEpipMap1 = aEpipModel.EpipMap1();
    const auto& anEpipMap2 = aEpipModel.EpipMap2();

    StdOut() << "Interpolator: " << aInterp->VNames() << ", Kernel Size: " << aInterp->SzKernel() << std::endl;

    auto aName1 = LastPrefix(FileOfPath(mNameIm1,false));
    auto aName2 = LastPrefix(FileOfPath(mNameIm2,false));
    auto anEpip1Name = aOutDir + replaceFirstOccurrence(replaceFirstOccurrence(mOutNamePat,"%1",aName1),"%2",aName2);
    auto anEpip2Name = aOutDir + replaceFirstOccurrence(replaceFirstOccurrence(mOutNamePat,"%1",aName2),"%2",aName1);
    auto aRPC1Name = anEpip1Name + ".xml";
    auto aRPC2Name = anEpip2Name + ".xml";

    // Resample Img1
    const auto* aIm1 = ReadIm2DGen(mNameIm1);
    StdOut() << "[EpipolarResample] Resampling image 1 "<< anEpipMap1.EpipFrame() << " Sz: " << anEpipMap1.EpipImSz() << std::endl;
    auto aIm1Rectif = aIm1->AllocReSampleGen(*aInterp, anEpipMap1, cTplBox(anEpipMap1.EpipImSz()));
    // Resample Img2
    const auto* aIm2 = ReadIm2DGen(mNameIm2);
    StdOut() << "[EpipolarResample] Resampling image 2 "<< anEpipMap2.EpipFrame() << " Sz: " << anEpipMap2.EpipImSz() << std::endl;
    auto aIm2Rectif = aIm2->AllocReSampleGen(*aInterp, anEpipMap2, cTplBox(anEpipMap2.EpipImSz()));

    StdOut() << "Image_1: " << anEpip1Name << std::endl;
    aIm1Rectif->ToFile(anEpip1Name);
    StdOut() << "RPC_1: " << aRPC1Name << std::endl;
    auto aResampSI1 = aSI1->GenerateSensorRPC( &anEpipMap1, nullptr, false, anEpip1Name);
    aResampSI1->ToFile(aRPC1Name);

    StdOut() << "Image2: " << anEpip2Name << std::endl;
    aIm2Rectif->ToFile(anEpip2Name);
    StdOut() << "RPC_2: " << aRPC2Name << std::endl;
    auto aResampSI2 = aSI2->GenerateSensorRPC(&anEpipMap2, nullptr, false, anEpip2Name );
    aResampSI2->ToFile(aRPC2Name);


    delete aResampSI1;
    delete aResampSI2;
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


cCollecSpecArg2007 & cAppli_EpipResampling::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    return anArgOpt
           << AOpt2007(mDegree,"Degree","Poly degree",{eTA2007::HDV})
           << AOpt2007(mDegreeInv,"DegreeInv","Inv Poly degree",{eTA2007::HDV})
           << AOpt2007(mNbByXY,"XYSteps","Nb XY steps",{eTA2007::HDV})
           << AOpt2007(mNbByZ,"ZSteps","Nb Z steps",{eTA2007::HDV})
           << AOpt2007(mMargin,"Margin","Output image margin",{eTA2007::HDV})
           << AOpt2007(mFrame,"FrameAlgo","Output image height algo",{eTA2007::HDV})
           << AOpt2007(mOutDir,"OutDir","Output directory (Default: VISU/" + Specs().Name()+")")
           << AOpt2007(mOutNamePat,"OutName","Output name pattern", {eTA2007::HDV})
           << AOpt2007(mInterpol,"Interpol","Interpolator", Append(cSpecOneArg2007::tAllSemPL{eTA2007::HDV},InterpolArgSem()))
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
        {eApF::ImProc},
        {eApDT::Orient,eApDT::Image},
        {eApDT::Orient,eApDT::Image},
        __FILE__
        );



/* ==================================================== */
/*                  TESTS                               */
/* ==================================================== */


class cAppli_EpipTest : public cMMVII_Appli
{
public :

    cAppli_EpipTest(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli &);
    int Exe() override;
    cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
    cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;

private :
    cPhotogrammetricProject  mPhProj;
    std::string  mNameIm1;
    cPt3dr mP;
};

cAppli_EpipTest::cAppli_EpipTest (
    const std::vector<std::string> &  aVArgs,
    const cSpecMMVII_Appli & aSpec
    )
    : cMMVII_Appli  (aVArgs,aSpec)
    , mPhProj       (*this)
{
}


cCollecSpecArg2007 & cAppli_EpipTest::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return anArgObl
           << Arg2007(mNameIm1,"Image",{eTA2007::FileImage})
           << mPhProj.DPOrient().ArgDirInMand()
           << Arg2007(mP,"Image Point",{})
        ;
}


cCollecSpecArg2007 & cAppli_EpipTest::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{

    return anArgOpt
        ;
}


int cAppli_EpipTest::Exe()
{
    mPhProj.FinishInit();


    const cSensorImage *  aSI1 =  mPhProj.ReadSensor(FileOfPath(mNameIm1,false /* Ok Not Exist*/),true/*DelAuto*/,false /* Not SVP*/);
    if (aSI1->HasIntervalZ()  && mP.z() == -1) {
        mP = cPt3dr(mP.x(), mP.y(), (aSI1->GetIntervalZ().x() + aSI1->GetIntervalZ().y()) / 2);
    }

    auto PG  = aSI1->ImageAndZ2Ground(mP);
    StdOut() << std::setprecision(10) << "PI: " << mP << " -> PG: " << PG << " -> PI:" << aSI1->Ground2ImageAndZ(PG)  << std::endl;
    return 0;

}


tMMVII_UnikPApli Alloc_EpipTest(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
    return tMMVII_UnikPApli(new cAppli_EpipTest(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_EpipTest
    (
        "EpipTest",
        Alloc_EpipTest,
        "Epipolar geometry tests (temporary for joe)",
        {eApF::Test},
        {eApDT::Orient,eApDT::Image},
        {eApDT::Orient,eApDT::Image},
        __FILE__
        );



}; // MMVII

