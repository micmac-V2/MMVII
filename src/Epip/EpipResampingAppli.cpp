#include "cMMVII_Appli.h"
#include "MMVII_Sensor.h"
#include "cEpipolarRectification.h"
#include "MMVII_Interpolators.h"
#include "MMVII_CodeTiming.h"
#include <vector>
#include <optional>
#include <cmath>

/**
   \file EpipGeom.cpp


 */


namespace MMVII
{

// Shared default pattern : EpipResampling reproduces it from the model's image names.
static constexpr const char* TheDefaultOutNamePat = "Epip_%1_%2.tif";

class cAppli_EpipRectification : public cMMVII_Appli
{
public :

    cAppli_EpipRectification(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli &);
    int Exe() override;
    cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
    cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;

private :
    void Resample(const std::string& aMasterName,
                  const std::string& aSlaveName,
                  const cSensorImage* aSI,
                  const cInterpolator1D* aInterp,
                  const cEpipolarMapping& anEpipMap
                  );

    cPhotogrammetricProject  mPhProj;
    std::string  mNameIm1;
    std::string  mNameIm2;
    int mDegree = 5;
    int mDegreeInv = mDegree + 4;
    int mNbByZ = 5;
    cPt2dr mZIntv;
    tREAL8 mTiePMaxRes = 2.0;
    tREAL8 mZMargin = 0.10;
    tREAL8 mTiePMinNbRatio = 0.04;
    int mTiePMinNbFloor = 25;
    tREAL8 mMaxResid = 0.1;
    std::string mOutDir;
    std::string mOutNamePat = TheDefaultOutNamePat;
    std::vector<std::string> mInterpol = {"Cubic","-0.5"};
    eEpipFrm mFrame = eEpipFrm::eIntersect;
    bool mSaveModel = false;
    bool mNoRPC = false;
    bool mNoImage = false;
};

cAppli_EpipRectification::cAppli_EpipRectification (
    const std::vector<std::string> &  aVArgs,
    const cSpecMMVII_Appli & aSpec
    )
    : cMMVII_Appli  (aVArgs,aSpec)
    , mPhProj       (*this)
{
}


// TODOCM: Gestion grosses images : daller ... Cache pour bout d'images ?
// TODOCM: Test d'epipolarisabilite ...
// TODOCM: Make output image name generation accessible for other apps
// TODOCM: Make sure image  file extension is present ? (.tif)

void cAppli_EpipRectification::Resample(const std::string& aMasterName,
                                     const std::string& aSlaveName,
                                     const cSensorImage* aSI,
                                     const cInterpolator1D* aInterp,
                                     const cEpipolarMapping& anEpipMap
                                     )
{
    auto aName1 = LastPrefix(FileOfPath(aMasterName,false));
    auto aName2 = LastPrefix(FileOfPath(aSlaveName,false));

    auto aBaseName = replaceFirstOccurrence(replaceFirstOccurrence(mOutNamePat,"%1",aName1),"%2",aName2);
    auto anEpipName = mOutDir + aBaseName;

    if (! mNoImage)
    {
        const auto* aIm = ReadIm2DGen(aMasterName);
        auto aIm1Rectif = aIm->AllocReSampleGen(*aInterp, anEpipMap, cTplBox(anEpipMap.EpipImSz()));
        StdOut() << "Name: " << anEpipName << std::endl;
        StdOut() << "Size: " << anEpipMap.EpipImSz() << std::endl;
        aIm1Rectif->ToFile(anEpipName);
        delete aIm1Rectif;
        delete aIm;
    }

    if (! mNoRPC)
    {
        auto aRPCName = mOutDir + "RPC_" + aBaseName + ".xml";
        StdOut() << "RPC : " << aRPCName << std::endl;
        auto aResampSI = aSI->GenerateSensorRPC( &anEpipMap, nullptr, false, anEpipName, anEpipMap.ZInterval());
        aResampSI->ToFile(aRPCName);
        delete aResampSI;
    }
}


int cAppli_EpipRectification::Exe()
{
    mPhProj.FinishInit();

    MMVII_INTERNAL_ASSERT_User(! mOutNamePat.empty(), eTyUEr::eUnClassedError, "OutName must not be empty");

    if (! IsInit(&mDegreeInv))
        mDegreeInv = mDegree + 4;

    if (! IsInit(&mOutDir))
    {
        mOutDir = mPhProj.DirVisuAppli();;
    }
    if (! mOutDir.empty())
    {
        mOutDir += "/";
    }
    CreateDirectories(mOutDir);
    const cInterpolator1D* aInterp = cDiffInterpolator1D::AllocFromNames(mInterpol);


    const cSensorImage *  aSI1 =  mPhProj.ReadSensor(FileOfPath(mNameIm1,false /* Ok Not Exist*/),true/*DelAuto*/,false /* Not SVP*/);
    const cSensorImage *  aSI2 =  mPhProj.ReadSensor(FileOfPath(mNameIm2,false /* Ok Not Exist*/),true/*DelAuto*/,false /* Not SVP*/);
    // Missing Z interval is an error unless ZIntv is given (checked once per sensor
    // in cEpipolarRectification); ZIntv overriding an existing interval only warns.
    const std::optional<cPt2dr> aZIntvArg = IsInit(&mZIntv) ? std::optional<cPt2dr>(mZIntv) : std::nullopt;

    // Create early to have possible error reported before doing long computations
    auto aDIm1 = cDataFileIm2D::Create(mNameIm1,eForceGray::No);
    auto aDIm2 = cDataFileIm2D::Create(mNameIm2,eForceGray::No);

    StdOut() << Color::sub_title << "*** Inputs" << Color::end << std::endl;
    StdOut() <<  "Image_1: " <<  mNameIm1;
    StdOut() << " " << aDIm1.Sz() << " " << ToStr(aDIm1.Type()) << " " << aDIm1.NbChannel() << " chan" << std::endl;
    StdOut() <<  "Image_2: " <<  mNameIm2;
    StdOut() << " " << aDIm2.Sz() << " "  << ToStr(aDIm2.Type()) << " " << aDIm2.NbChannel() << " chan" << std::endl;

    StdOut() << "Degree: " << mDegree << ", DegreeInv: " << mDegreeInv << std::endl;
    StdOut() << "NbByZ: " << mNbByZ << std::endl;
    StdOut() << "Frame: " << ToStr(mFrame) << std::endl;
    StdOut() << "Interpolator: " << aInterp->VNames() << ", Kernel Size: " << aInterp->SzKernel() << std::endl;

    StdOut() << Color::sub_title << "*** Rectification" << Color::end << std::endl;
    auto aParams = cEpipolarRectification::cParams{mDegree,mDegreeInv,mNbByZ,mFrame};
    aParams.mZIntv = aZIntvArg;
    aParams.mTiePMaxRes = mTiePMaxRes;
    aParams.mZMargin = mZMargin;
    aParams.mTiePMinNbRatio = mTiePMinNbRatio;
    aParams.mTiePMinNbFloor = mTiePMinNbFloor;
    if (mPhProj.DPTieP().DirInIsInit())
    {
        cSetHomogCpleIm aSetH;
        bool aHasHom = mPhProj.GenReadHomol(aSetH, mNameIm1, mNameIm2);
        MMVII_INTERNAL_ASSERT_User(aHasHom && (aSetH.NbH() > 0), eTyUEr::eOpenFile,
            "No tie points found between the two images in the TieP directory");
        aParams.mHomolPts = aSetH;
    }
    auto aRectifier = cEpipolarRectification(*aSI1, *aSI2, aParams);
    auto aEpipModel = aRectifier.Compute();

    StdOut() << "Nb Pairs 1->2 : " << aRectifier.NbPairs12() << std::endl;
    StdOut() << "Nb Pairs 2->1 : " << aRectifier.NbPairs21() << std::endl;

    for (const auto& [aName,aMap] : {std::make_pair("Image_1",&aEpipModel.EpipMap1()), std::make_pair("Image_2",&aEpipModel.EpipMap2())})
    {
        StdOut() << "Grid " << aName << " : step=" << aMap->GridStep() << "px, "
                 << aMap->NbStepX() << "*" << aMap->NbStepY() << "=" << (aMap->NbStepX()*aMap->NbStepY()) << " cells" << std::endl;
    }

    // Independent (held-out) residual check, complementing the train-biased variance above.
    const tREAL8 aV1V2ResidIndep = std::sqrt(aRectifier.V1V2VarIndep());
    const tREAL8 aW1ResidIndep   = std::sqrt(aRectifier.W1VarIndep());
    const tREAL8 aW2ResidIndep   = std::sqrt(aRectifier.W2VarIndep());
    StdOut() << "V1,V2 errors sigma (indep, px) : " << Color::info << aV1V2ResidIndep << Color::end << std::endl;
    StdOut() << "W1 errors sigma (indep, px) : " << Color::info << aW1ResidIndep << Color::end << std::endl;
    StdOut() << "W2 errors sigma (indep, px) : " << Color::info << aW2ResidIndep << Color::end << std::endl;
    if (aV1V2ResidIndep > mMaxResid)
    {
        MMVII_UserError(eTyUEr::eUnClassedError,
            "Independent V1/V2 residual too high (" + ToStr(aV1V2ResidIndep) + " > " + ToStr(mMaxResid) + ")");
    }
    if (aW1ResidIndep > mMaxResid)
    {
        MMVII_UserError(eTyUEr::eUnClassedError,
            "Independent W1 residual too high (" + ToStr(aW1ResidIndep) + " > " + ToStr(mMaxResid) + ")");
    }
    if (aW2ResidIndep > mMaxResid)
    {
        MMVII_UserError(eTyUEr::eUnClassedError,
            "Independent W2 residual too high (" + ToStr(aW2ResidIndep) + " > " + ToStr(mMaxResid) + ")");
    }


    const auto& anEpipMap1 = aEpipModel.EpipMap1();
    const auto& anEpipMap2 = aEpipModel.EpipMap2();

    auto aName1 = LastPrefix(FileOfPath(mNameIm1,false));
    auto aName2 = LastPrefix(FileOfPath(mNameIm2,false));

    if (mSaveModel)
    {
        // One file per image, named like the resampled image on that side.
        StdOut() << Color::sub_title << "*** Model" << Color::end << std::endl;
        const std::string & aOriName = mPhProj.DPOrient().DirIn();

        auto aBaseName1 = replaceFirstOccurrence(replaceFirstOccurrence(mOutNamePat,"%1",aName1),"%2",aName2);
        auto aModelName1 = mOutDir + LastPrefix(aBaseName1) + ".EpipModel." + GlobTaggedNameDefSerial();
        cEpipSideModel(static_cast<const cEpipPolyMapping&>(anEpipMap1), aOriName, mNameIm1, mNameIm2).ToFile(aModelName1);
        StdOut() << "Image_1: " << aModelName1 << std::endl;

        auto aBaseName2 = replaceFirstOccurrence(replaceFirstOccurrence(mOutNamePat,"%1",aName2),"%2",aName1);
        auto aModelName2 = mOutDir + LastPrefix(aBaseName2) + ".EpipModel." + GlobTaggedNameDefSerial();
        cEpipSideModel(static_cast<const cEpipPolyMapping&>(anEpipMap2), aOriName, mNameIm2, mNameIm1).ToFile(aModelName2);
        StdOut() << "Image_2: " << aModelName2 << std::endl;
    }

    if ((! mNoImage) || (! mNoRPC))
    {
        StdOut() << Color::sub_title << "*** Resampling" << Color::end << std::endl;

        // Resample Img1
        StdOut() << Color::title << "* Image 1" << Color::end << std::endl;
        Resample(mNameIm1,mNameIm2,aSI1,aInterp,anEpipMap1);

        // Resample Img2
        StdOut() << Color::title << "* Image 2" << Color::end << std::endl;
        Resample(mNameIm2,mNameIm1,aSI2,aInterp,anEpipMap2);
    }


    delete aInterp;
    return EXIT_SUCCESS;
}


cCollecSpecArg2007 & cAppli_EpipRectification::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return anArgObl
          << Arg2007(mNameIm1,"name first image",{eTA2007::FileImage})
          << Arg2007(mNameIm2,"name second image",{eTA2007::FileImage})
          << mPhProj.DPOrient().ArgDirInMand()
        ;
}


cCollecSpecArg2007 & cAppli_EpipRectification::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    return anArgOpt
           << cHeaderSectionArg("Rectification")
           << AOpt2007(mDegree,"Degree","Poly degree",{eTA2007::HDV})
           << AOpt2007(mDegreeInv,"DegreeInv","Inv Poly degree",{eTA2007::HDV})
           << AOpt2007(mMaxResid,"MaxResid","Max independent (held-out) V1/V2, W1 and W2 residual (px), error above",{eTA2007::HDV})
           << cHeaderSectionArg("Z interval")
           << AOpt2007(mNbByZ,"ZSteps","Nb Z steps",{eTA2007::HDV})
           << AOpt2007(mZIntv,"ZIntv","Z interval [Zmin,Zmax], overrides sensor's own and any TieP-derived one (mandatory when sensor has none, e.g. no RPC, and TieP is not given either)")
           << cHeaderSectionArg("Tie Points")
           << mPhProj.DPTieP().ArgDirInOpt("TieP","Tie points to infer Z interval from (alternative to ZIntv, overrides sensor's own)")
           << AOpt2007(mZMargin,"ZMargin","Relative margin added around the raw [Zmin,Zmax] envelope inferred from TieP",{eTA2007::HDV})
           << AOpt2007(mTiePMaxRes,"TiePMaxRes","Max triangulation residual (px) for a tie point kept when inferring Z from TieP",{eTA2007::HDV})
           << AOpt2007(mTiePMinNbRatio,"TiePMinNbRatio","Min kept tie points = max(TiePMinNbFloor,ratio*sqrt(W*H))",{eTA2007::HDV})
           << AOpt2007(mTiePMinNbFloor,"TiePMinNbFloor","Absolute floor for the min kept tie point count",{eTA2007::HDV})
           << cHeaderSectionArg("Resampling")
           << AOpt2007(mInterpol,"Interpol","Interpolator", Append(cSpecOneArg2007::tAllSemPL{eTA2007::HDV},InterpolArgSem()))
           << AOpt2007(mFrame,"FrameAlgo","Output image height algo",{eTA2007::HDV})
           << cHeaderSectionArg("Output")
           << AOpt2007(mOutDir,"OutDir","Output directory (Default: VISU/" + Specs().Name()+")")
           << AOpt2007(mOutNamePat,"OutName","Output name pattern for images and other output files (must not be empty)", {eTA2007::HDV})
           << AOpt2007(mSaveModel,"SaveModel","Serialize the computed model (name derived from OutName)",{eTA2007::HDV})
           << AOpt2007(mNoRPC,"NoRPC","Don't (re)compute the resampled RPC",{eTA2007::HDV})
           << AOpt2007(mNoImage,"NoImage","Don't produce the resampled TIF images",{eTA2007::HDV})
        ;
}



/* ==================================================== */

tMMVII_UnikPApli Alloc_EpipRectification(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
    return tMMVII_UnikPApli(new cAppli_EpipRectification(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_EpipRectification
    (
        "EpipRectification",
        Alloc_EpipRectification,
        "Compute the epipolar geometry of two images, and optionally resample them",
        {eApF::ImProc},
        {eApDT::Orient,eApDT::Image},
        {eApDT::Orient,eApDT::Image},
        __FILE__
        );



/* ==================================================== */
/*         EpipResampling : resampling a posteriori      */
/* ==================================================== */

// Adds a crop+scale on top of a cEpipolarMapping, at resampling time only (never serialized).
class cEpipCropScaleMapping : public cDataInvertibleMapping<tREAL8,2>
{
public:
    cEpipCropScaleMapping(const cEpipolarMapping& aMap, cPt2dr aCropP0, tREAL8 aScale)
        : mMap(aMap), mCropP0(aCropP0), mScale(aScale) {}

    cPt2dr Value(const cPt2dr& aPt) const override
    {
        return (mMap.Value(aPt) - mCropP0) * mScale;
    }
    cPt2dr Inverse(const cPt2dr& aPt) const override
    {
        return mMap.Inverse(aPt / mScale + mCropP0);
    }

private:
    const cEpipolarMapping& mMap;
    cPt2dr mCropP0;
    tREAL8 mScale;
};


class cAppli_EpipResampling : public cMMVII_Appli
{
public :

    cAppli_EpipResampling(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli &);
    int Exe() override;
    cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
    cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;

private :
    cPhotogrammetricProject  mPhProj;
    std::string  mNameIm;
    std::string  mNameModel;
    std::vector<std::string> mInterpol = {"Cubic","-0.5"};
    std::string  mOutDir;
    std::string  mOutName;
    cPt2di       mCropP0{0,0};
    cPt2di       mCropP1{0,0};
    tREAL8       mScale = 1.0;
    bool         mNoRPC = false;
    bool         mNoImage = false;
};

cAppli_EpipResampling::cAppli_EpipResampling (
    const std::vector<std::string> &  aVArgs,
    const cSpecMMVII_Appli & aSpec
    )
    : cMMVII_Appli  (aVArgs,aSpec)
    , mPhProj       (*this)
{
}

int cAppli_EpipResampling::Exe()
{
    mPhProj.FinishInit();
    auto aSideModel = cEpipSideModel::FromFile(mNameModel);

    if (! IsInit(&mOutDir))
        mOutDir = mPhProj.DirVisuAppli();
    if (! mOutDir.empty())
        mOutDir += "/";
    CreateDirectories(mOutDir);

    const cEpipolarMapping & anEpipMap = aSideModel.Map();

    MMVII_INTERNAL_ASSERT_User(IsInit(&mCropP0) == IsInit(&mCropP1), eTyUEr::eUnClassedError,
        "CropP0 and CropP1 must be given together");
    const cPt2dr aCropP0 = IsInit(&mCropP0) ? ToR(mCropP0) : cPt2dr(0,0);
    const cPt2dr aCropP1 = IsInit(&mCropP1) ? ToR(mCropP1) : ToR(anEpipMap.EpipImSz());

    cEpipCropScaleMapping aResampMap(anEpipMap, aCropP0, mScale);
    const cPt2di aOutSz( round_ni((aCropP1.x() - aCropP0.x()) * mScale),
                         round_ni((aCropP1.y() - aCropP0.y()) * mScale) );

    const cInterpolator1D* aInterp = cDiffInterpolator1D::AllocFromNames(mInterpol);

    // Reproduces EpipRectification's own default naming from the model's image names.
    const std::string aDefaultOutName = replaceFirstOccurrence(replaceFirstOccurrence(
            std::string(TheDefaultOutNamePat),
            "%1",LastPrefix(FileOfPath(aSideModel.MasterImName(),false))),
            "%2",LastPrefix(FileOfPath(aSideModel.SlaveImName(),false)));
    const std::string aOutBaseName = IsInit(&mOutName) ? mOutName : aDefaultOutName;
    auto anOutName = mOutDir + aOutBaseName;

    StdOut() << "Name: " << anOutName << std::endl;

    if (! mNoImage)
    {
        StdOut() << "Size: " << aOutSz << std::endl;
        const auto* aIm = ReadIm2DGen(mNameIm);
        auto aImRectif = aIm->AllocReSampleGen(*aInterp, aResampMap, cTplBox(aOutSz));
        aImRectif->ToFile(anOutName);
        delete aImRectif;
        delete aIm;
    }

    cSensorImage* aResampSI = nullptr;
    if (! mNoRPC)
    {
        // Two-pass pattern (cf. cAppliMeshImageDevlp::Exe()) : SetDirIn needs a 2nd FinishInit().
        mPhProj.DPOrient().SetDirIn(aSideModel.OriName());
        mPhProj.FinishInit();
        const cSensorImage * aSI = mPhProj.ReadSensor(FileOfPath(mNameIm,false /* Ok Not Exist*/),true/*DelAuto*/,false /* Not SVP*/);

        auto aRPCName = mOutDir + "RPC_" + aOutBaseName + ".xml";
        StdOut() << "RPC : " << aRPCName << std::endl;
        aResampSI = aSI->GenerateSensorRPC(&aResampMap, nullptr, false, anOutName, anEpipMap.ZInterval());
        aResampSI->ToFile(aRPCName);
    }

    delete aResampSI;
    delete aInterp;
    return EXIT_SUCCESS;
}

cCollecSpecArg2007 & cAppli_EpipResampling::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return anArgObl
          << Arg2007(mNameIm,"name of the image to resample",{eTA2007::FileImage})
          << Arg2007(mNameModel,"epipolar model file for that image, as saved by EpipRectification's SaveModel=true",{eTA2007::FileTagged})
        ;
}

cCollecSpecArg2007 & cAppli_EpipResampling::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    return anArgOpt
           << cHeaderSectionArg("Resampling")
           << AOpt2007(mInterpol,"Interpol","Interpolator", Append(cSpecOneArg2007::tAllSemPL{eTA2007::HDV},InterpolArgSem()))
           << AOpt2007(mCropP0,"CropP0","Crop : lower-left corner in the model's epipolar output coordinates (with CropP1 ; default : full frame)")
           << AOpt2007(mCropP1,"CropP1","Crop : upper-right corner in the model's epipolar output coordinates (with CropP0 ; default : full frame)")
           << AOpt2007(mScale,"Scale","Resample at a different pixel density over the (possibly cropped) region ; >1 zooms in, <1 zooms out",{eTA2007::HDV})
           << AOpt2007(mNoRPC,"NoRPC","Don't (re)compute the resampled RPC ; also skips reopening the sensor",{eTA2007::HDV})
           << AOpt2007(mNoImage,"NoImage","Don't produce the resampled TIF image (e.g. to only regenerate its RPC)",{eTA2007::HDV})
           << cHeaderSectionArg("Output")
           << AOpt2007(mOutDir,"OutDir","Output directory (Default: VISU/" + Specs().Name()+")")
           << AOpt2007(mOutName,"OutName","Output image name (Default: same pattern/name EpipRectification would have used)")
        ;
}


tMMVII_UnikPApli Alloc_EpipResampling(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
    return tMMVII_UnikPApli(new cAppli_EpipResampling(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_EpipResampling
    (
        "EpipResampling",
        Alloc_EpipResampling,
        "Resample one image in epipolar geometry from a model saved by EpipRectification (SaveModel=true), with optional crop/scale",
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

