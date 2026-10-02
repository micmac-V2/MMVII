#include "cMMVII_Appli.h"
#include "MMVII_Sensor.h"
#include "cEpipolarRectification.h"
#include "MMVII_Interpolators.h"
#include "MMVII_CodeTiming.h"
#include "../Sensors/cExternalSensor.h"
#include <vector>
#include <array>
#include <optional>
#include <cmath>

/**
   \file EpipGeom.cpp


 */


namespace MMVII
{

// Shared default pattern, used by both commands from the image names.
static constexpr const char* TheDefaultOutNamePat = "Epip_%1_%2.tif";

// Output base name of an image : %1 = this image, %2 = the other image of the pair.
static std::string EpipOutBaseName(const std::string & aPattern,const std::string & aNameIm,const std::string & aNameOther)
{
    return replaceFirstOccurrence(replaceFirstOccurrence(aPattern,
            "%1",LastPrefix(FileOfPath(aNameIm,false))),
            "%2",LastPrefix(FileOfPath(aNameOther,false)));
}
// Info file of a pair of crops, named from the output name of the master image
static std::string EpipInfoName(const std::string & aOutDir,const std::string & aMasterBaseName)
{
    return aOutDir + LastPrefix(aMasterBaseName) + ".Info." + GlobTaggedNameDefSerial();
}

// Crop options of image 1 and 2 : the crop is given on image anOpt.mMaster, the other one is derived from the Z interval.
// Writes the info file of the pair of crops (images aOutBaseNames[0..1]) when aWriteInfo, i.e. when the images are resampled ;
// without a crop, the crops are the whole frames. Without crop nor info file there is nothing to derive.
static std::array<cEpipCropMaskOpts,2> EpipPairCropOpts(const cEpipCropMaskOpts & anOpt,
        const cEpipolarMapping & aMap1, const cEpipolarMapping & aMap2,
        const cSensorImage & aSI1, const cSensorImage & aSI2,
        const std::string & aOutDir, const std::array<std::string,2> & aOutBaseNames, bool aWriteInfo)
{
    std::array<cEpipCropMaskOpts,2> aRes{anOpt,anOpt};
    if ((! anOpt.mHasCrop) && (! aWriteInfo))
        return aRes;
    MMVII_INTERNAL_ASSERT_User((anOpt.mMaster==1) || (anOpt.mMaster==2), eTyUEr::eUnClassedError, "Master must be 1 or 2");

    const int aKM = anOpt.mMaster - 1;   // master and slave indices in 0..1
    const cEpipolarMapping * aMaps[2] = {&aMap1,&aMap2};
    const cSensorImage * aSIs[2] = {&aSI1,&aSI2};
    const cPt2dr aZIntv = EpipPairZInterval(aMap1,aMap2);
    const cPt2di aMasterP0 = anOpt.mHasCrop ? anOpt.mCropP0 : cPt2di(0,0);
    const cPt2di aMasterP1 = anOpt.mHasCrop ? anOpt.mCropP1 : aMaps[aKM]->EpipImSz();
    cEpipSlaveCrop aSlave = EpipSlaveCrop(*aMaps[aKM],*aMaps[1-aKM],*aSIs[aKM],*aSIs[1-aKM],aMasterP0,aMasterP1,aZIntv);
    if (anOpt.mHasCrop)
    {
        aRes[1-aKM].mCropP0 = aSlave.mP0;
        aRes[1-aKM].mCropP1 = aSlave.mP1;
        StdOut() << "Crop image " << 2-aKM << " (derived) : " << aSlave.mP0 << " " << aSlave.mP1 << std::endl;
    }
    else
    {
        aSlave.mP0 = cPt2di(0,0);   // no crop : the whole slave frame, the disparity range is that of the whole master frame
        aSlave.mP1 = aMaps[1-aKM]->EpipImSz();
    }

    StdOut() << "Z interval : " << aZIntv << std::endl;
    StdOut() << "Disparity range (slave - master, epipolar coordinates) : " << aSlave.mDispRange << std::endl;
    // Parallax due to Z alone (the range also varies with the position) : under 1 px, no depth information
    // (nearly parallel views, or a Z interval too small) ; over half the master frame, almost surely a wrong Z interval
    StdOut() << "Mean parallax over the Z interval : " << aSlave.mMeanParallax << " px" << std::endl;
    if (aSlave.mMeanParallax < 1.0)
        MMVII_USER_WARNING("The parallax over the Z interval is under 1 px (" + ToStr(aSlave.mMeanParallax)
                           + ") : no depth information (nearly parallel views, or a Z interval too small)");
    if (aSlave.mMeanParallax > 0.5 * aMaps[aKM]->EpipImSz().x())
        MMVII_USER_WARNING("The parallax over the Z interval (" + ToStr(aSlave.mMeanParallax)
                           + " px) is wider than half of the master epipolar image : check the Z interval");

    if (aWriteInfo)
    {
        const cEpipCropInfo anInfo = MakeEpipCropInfo(aMasterP0,aMasterP1,aSlave,aZIntv,aMaps[aKM]->EpipImSz(),aMaps[1-aKM]->EpipImSz(),
                                                      aOutBaseNames[aKM],aOutBaseNames[1-aKM]);
        const std::string anInfoName = EpipInfoName(aOutDir,aOutBaseNames[aKM]);
        anInfo.ToFile(anInfoName);
        StdOut() << "Info : " << anInfoName << std::endl;
    }
    return aRes;
}


// Resampling map of the crop expressed in a window of the source image (origin of the window : aWinP0)
class cEpipWindowMapping : public cDataInvertibleMapping<tREAL8,2>
{
public:
    cEpipWindowMapping(const cDataInvertibleMapping<tREAL8,2> & aMap, const cPt2di & aWinP0)
        : mMap(aMap), mWinP0(ToR(aWinP0)) {}

    cPt2dr Value(const cPt2dr& aPt) const override { return mMap.Value(aPt + mWinP0); }
    cPt2dr Inverse(const cPt2dr& aPt) const override { return mMap.Inverse(aPt) - mWinP0; }

private:
    const cDataInvertibleMapping<tREAL8,2> & mMap;
    cPt2dr mWinP0;
};

// Window of the source (clipped to it) read to resample a crop : footprint of the output pixels
// (sampled on a grid) widened by the interpolation kernel. False if the crop is outside the source.
static bool SourceWindow(const cDataInvertibleMapping<tREAL8,2> & aMap, const cPt2di & aOutSz,
                         const cPt2di & aSzSrc, tREAL8 aSzKernel, cBox2di & aWin)
{
    const int aNbSamp = 33;
    cPt2dr aMin(1e30,1e30), aMax(-1e30,-1e30);
    for (int aKx=0 ; aKx<aNbSamp ; aKx++)
        for (int aKy=0 ; aKy<aNbSamp ; aKy++)
        {
            const cPt2dr aPIn = aMap.Inverse(cPt2dr((aOutSz.x()-1)*aKx/(aNbSamp-1.0),(aOutSz.y()-1)*aKy/(aNbSamp-1.0)));
            aMin = cPt2dr(std::min(aMin.x(),aPIn.x()),std::min(aMin.y(),aPIn.y()));
            aMax = cPt2dr(std::max(aMax.x(),aPIn.x()),std::max(aMax.y(),aPIn.y()));
        }
    const int aMarg = round_up(aSzKernel) + 2;   // kernel + safety for the grid sampling
    const cPt2di aP0(std::max(0,round_down(aMin.x())-aMarg), std::max(0,round_down(aMin.y())-aMarg));
    const cPt2di aP1(std::min(aSzSrc.x(),round_down(aMax.x())+1+aMarg), std::min(aSzSrc.y(),round_down(aMax.y())+1+aMarg));
    if ((aP1.x()<=aP0.x()) || (aP1.y()<=aP0.y()))
        return false;
    aWin = cBox2di(aP0,aP1);
    return true;
}

void ResampleEpipImage(const cEpipCropMaskOpts & anOpt,
                              const std::string & aNameIm,
                              const cEpipolarMapping & anEpipMap,
                              const cSensorImage * aSI,
                              const std::string & aRPCFile,
                              const cInterpolator1D & aInterp,
                              bool aNoImage,
                              const std::string & aOutDir,
                              const std::string & aOutBaseName)
{
    const cPt2dr aCropP0 = anOpt.mHasCrop ? ToR(anOpt.mCropP0) : cPt2dr(0,0);
    const cPt2dr aCropP1 = anOpt.mHasCrop ? ToR(anOpt.mCropP1) : ToR(anEpipMap.EpipImSz());
    cEpipCropMapping aResampMap(anEpipMap, aCropP0);
    const cPt2di aOutSz( round_ni(aCropP1.x() - aCropP0.x()), round_ni(aCropP1.y() - aCropP0.y()) );
    MMVII_INTERNAL_ASSERT_User((aOutSz.x() > 0) && (aOutSz.y() > 0), eTyUEr::eUnClassedError,
        "CropP1 must be strictly above/right of CropP0");

    const std::string anOutName = aOutDir + aOutBaseName;
    const cDataFileIm2D aDFSrc = cDataFileIm2D::Create(aNameIm,eForceGray::Yes);

    // A pixel is valid when it maps inside the source image, i.e. the criterion AllocReSampleGen uses to interpolate.
    if (anOpt.mMaskOn)
    {
        const std::string aMaskName = aOutDir + replaceFirstOccurrence(anOpt.mMaskName,"$1",LastPrefix(aOutBaseName));
        const cPt2di aSzSrc = aDFSrc.Sz();
        cIm2D<tU_INT1> aMaskIm(aOutSz);
        for (const auto & aPix : aMaskIm.DIm())
            aMaskIm.DIm().SetV(aPix, cRect2(cPt2di(0,0),aSzSrc).Inside(ToI(aResampMap.Inverse(ToR(aPix)))) ? 1 : 0);
        StdOut() << "Mask: " << aMaskName << std::endl;
        aMaskIm.DIm().ToFile(aMaskName,eTyNums::eTN_U_INT1,{"NBITS=1"});
    }

    if (! aNoImage)
    {
        StdOut() << "Name: " << anOutName << std::endl;
        StdOut() << "Size: " << aOutSz << std::endl;
        cDataGenUnTypedIm<2> * aImRectif = nullptr;
        cBox2di aWin = cBox2di::Empty();
        if (SourceWindow(aResampMap,aOutSz,aDFSrc.Sz(),aInterp.SzKernel(),aWin))
        {
            // Only the window of the source needed by the crop is read
            const auto* aIm = ReadIm2DGen(aNameIm,aWin);
            aImRectif = aIm->AllocReSampleGen(aInterp, cEpipWindowMapping(aResampMap,aWin.P0()), cTplBox(aOutSz));
            delete aIm;
        }
        else
        {
            // Crop outside the source : null image
            aImRectif = AllocImGen(aOutSz,aDFSrc.Type());
            for (const auto & aPix : *aImRectif)
                aImRectif->VD_SetV(aPix,0);
        }
        aImRectif->ToFile(anOutName);
        delete aImRectif;
    }

    if (aSI)
    {
        auto aRPCName = aOutDir + "RPC_" + aOutBaseName + ".xml";
        StdOut() << "RPC : " << aRPCName << std::endl;
        // RPC of the full frame (read, or fitted), then cropped : all crops of an image share the same polynomials.
        auto aResampSI = aRPCFile.empty()
            ? aSI->GenerateSensorRPC(&anEpipMap, nullptr, false, anOutName, anEpipMap.ZInterval())
            : ReadExternalSensor(aRPCFile, anOutName, false);
        if (anOpt.mHasCrop)
        {
            auto aCropSI = aResampSI->CropSensor(ToI(aCropP0), aOutSz);
            delete aResampSI;
            aResampSI = aCropSI;
        }
        aResampSI->ToFile(aRPCName);
        delete aResampSI;
    }
}

class cAppli_EpipRectification : public cMMVII_Appli
{
public :

    cAppli_EpipRectification(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli &);
    int Exe() override;
    cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
    cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;
    std::vector<cOneHelpSampleCmp> Samples() const override;

private :
    void Resample(const std::string& aMasterName,
                  const std::string& aSlaveName,
                  const cSensorImage* aSI,
                  const cInterpolator1D* aInterp,
                  const cEpipolarMapping& anEpipMap,
                  const cEpipCropMaskOpts& aCropMask
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
    cEpipCropMaskOpts mCropMask;
};

cAppli_EpipRectification::cAppli_EpipRectification (
    const std::vector<std::string> &  aVArgs,
    const cSpecMMVII_Appli & aSpec
    )
    : cMMVII_Appli  (aVArgs,aSpec)
    , mPhProj       (*this)
{
}


void cAppli_EpipRectification::Resample(const std::string& aMasterName,
                                     const std::string& aSlaveName,
                                     const cSensorImage* aSI,
                                     const cInterpolator1D* aInterp,
                                     const cEpipolarMapping& anEpipMap,
                                     const cEpipCropMaskOpts& aCropMask
                                     )
{
    auto aBaseName = EpipOutBaseName(mOutNamePat,aMasterName,aSlaveName);
    ResampleEpipImage(aCropMask, aMasterName, anEpipMap, mNoRPC ? nullptr : aSI, "", *aInterp, mNoImage, mOutDir, aBaseName);
}


int cAppli_EpipRectification::Exe()
{
    mPhProj.FinishInit();

    MMVII_INTERNAL_ASSERT_User(! mOutNamePat.empty(), eTyUEr::eUnClassedError, "OutName must not be empty");
    // Names of the written images need an extension (GDAL chooses the format from it)
    mOutNamePat = EpipNameWithExtension(mOutNamePat);
    mCropMask.mMaskName = EpipNameWithExtension(mCropMask.mMaskName);

    mCropMask.mMaskOn = mCropMask.mMask || IsInit(&mCropMask.mMaskName);   // no crop here, image 1 is the master

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

    if (mSaveModel)
    {
        // One file for the pair, named like the first image's resampled output.
        StdOut() << Color::sub_title << "*** Model" << Color::end << std::endl;
        auto aBaseName = EpipOutBaseName(mOutNamePat,mNameIm1,mNameIm2);
        auto aModelName = mOutDir + LastPrefix(aBaseName) + ".EpipModel." + GlobTaggedNameDefSerial();
        cEpipPairModel aModel(static_cast<const cEpipPolyMapping&>(anEpipMap1), static_cast<const cEpipPolyMapping&>(anEpipMap2),
                              mPhProj.DPOrient().DirIn(), mNameIm1, mNameIm2);
        if (! mNoRPC)   // the RPC of the full frames are written with the images, in the directory of the model
            aModel.SetRPCNames("RPC_" + EpipOutBaseName(mOutNamePat,mNameIm1,mNameIm2) + ".xml",
                               "RPC_" + EpipOutBaseName(mOutNamePat,mNameIm2,mNameIm1) + ".xml");
        aModel.ToFile(aModelName);
        StdOut() << "Model: " << aModelName << std::endl;
    }

    if ((! mNoImage) || (! mNoRPC))
    {
        StdOut() << Color::sub_title << "*** Resampling" << Color::end << std::endl;

        // The info file describes the resampled images : written only when they are
        const auto aCropOpts = EpipPairCropOpts(mCropMask,anEpipMap1,anEpipMap2,*aSI1,*aSI2,mOutDir,
                {EpipOutBaseName(mOutNamePat,mNameIm1,mNameIm2),EpipOutBaseName(mOutNamePat,mNameIm2,mNameIm1)},! mNoImage);

        // Resample Img1
        StdOut() << Color::title << "* Image 1" << Color::end << std::endl;
        Resample(mNameIm1,mNameIm2,aSI1,aInterp,anEpipMap1,aCropOpts[0]);

        // Resample Img2
        StdOut() << Color::title << "* Image 2" << Color::end << std::endl;
        Resample(mNameIm2,mNameIm1,aSI2,aInterp,anEpipMap2,aCropOpts[1]);
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
    anArgOpt
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
        << AOpt2007(mFrame,"FrameAlgo","Output image height algo",{eTA2007::HDV})
        << AOpt2007(mInterpol,"Interpol","Interpolator", Append(cSpecOneArg2007::tAllSemPL{eTA2007::HDV},InterpolArgSem()))
        << cHeaderSectionArg("Output")
        << AOpt2007(mOutDir,"OutDir","Output directory (Default: VISU/" + Specs().Name()+")")
        << AOpt2007(mOutNamePat,"OutName","Output name pattern for images and other output files (must not be empty ; .tif is added if there is no extension)", {eTA2007::HDV})
        << AOpt2007(mSaveModel,"SaveModel","Serialize the computed model (name derived from OutName)",{eTA2007::HDV})
        << AOpt2007(mNoRPC,"NoRPC","Don't write the RPC of the resampled images",{eTA2007::HDV})
        << AOpt2007(mNoImage,"NoImage","Don't produce the resampled TIF images",{eTA2007::HDV})
        << AOpt2007(mCropMask.mMask,"Mask","Also write a 1-bit mask of the valid pixels of the resampled image (1 = maps inside the source image)",{eTA2007::HDV})
        << AOpt2007(mCropMask.mMaskName,"MaskName","Mask name pattern, $1 = output image name without extension ; implies Mask",{eTA2007::HDV})
        ;
    return anArgOpt;
}



std::vector<cOneHelpSampleCmp> cAppli_EpipRectification::Samples() const
{
    return
    {
        cOneHelpSampleCmp::Header("Rectify a pair (RPC sensors) and save the model"),
        {"MMVII EpipRectification Im1.tif Im2.tif Ori SaveModel=true"},
        cOneHelpSampleCmp::Header("Central perspective cameras : Z interval given or inferred from tie points"),
        {"MMVII EpipRectification Im1.tif Im2.tif Ori ZIntv=[0,100]"},
        {"MMVII EpipRectification Im1.tif Im2.tif Ori TieP=Std"},
        cOneHelpSampleCmp::Header("With validity masks (to crop, use EpipResampling on the saved model)"),
        {"MMVII EpipRectification Im1.tif Im2.tif Ori SaveModel=true Mask=true"},
        cOneHelpSampleCmp::Header("Only compute and save the model, no resampling"),
        {"MMVII EpipRectification Im1.tif Im2.tif Ori SaveModel=true NoImage=true NoRPC=true"}
    };
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

class cAppli_EpipResampling : public cMMVII_Appli
{
public :

    cAppli_EpipResampling(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli &);
    int Exe() override;
    cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
    cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;
    std::vector<cOneHelpSampleCmp> Samples() const override;

private :
    int ExeTiles(const cEpipPairModel & aModel);   ///< Cut the pair in tiles, one command per tile in parallel

    cPhotogrammetricProject  mPhProj;
    std::string  mNameModel;
    cPt2di       mSzTiles{0,0};
    cPt2di       mSzOverL{0,0};
    std::vector<std::string> mInterpol = {"Cubic","-0.5"};
    std::string  mOutDir;
    std::string  mOutNamePat = TheDefaultOutNamePat;
    bool         mNoRPC = false;
    cEpipCropMaskOpts mCropMask;
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
    auto aModel = cEpipPairModel::FromFile(mNameModel);

    MMVII_INTERNAL_ASSERT_User(! mOutNamePat.empty(), eTyUEr::eUnClassedError, "OutName must not be empty");
    // Names of the written images need an extension (GDAL chooses the format from it)
    mOutNamePat = EpipNameWithExtension(mOutNamePat);
    mCropMask.mMaskName = EpipNameWithExtension(mCropMask.mMaskName);

    if (! IsInit(&mOutDir))
        mOutDir = mPhProj.DirVisuAppli();
    if (! mOutDir.empty())
        mOutDir += "/";
    CreateDirectories(mOutDir);

    MMVII_INTERNAL_ASSERT_User(IsInit(&mCropMask.mCropP0) == IsInit(&mCropMask.mCropP1), eTyUEr::eUnClassedError,
        "CropP0 and CropP1 must be given together");
    mCropMask.mHasCrop = IsInit(&mCropMask.mCropP0);
    mCropMask.mMaskOn = mCropMask.mMask || IsInit(&mCropMask.mMaskName);

    // A tile of a tiling is itself a call of this command (level > 0) : it never cuts again
    if (IsInit(&mSzTiles) && (LevelCall()==0))
        return ExeTiles(aModel);

    const cInterpolator1D* aInterp = cDiffInterpolator1D::AllocFromNames(mInterpol);

    // Sensors are needed for the RPC, and for the info file (slave crop and disparity range from the Z interval)
    std::string aStoredNames[2], aNameIms[2];
    const cSensorImage* aSIs[2] = {nullptr,nullptr};
    // Two-pass pattern (cf. cAppliMeshImageDevlp::Exe()) : SetDirIn needs a 2nd FinishInit().
    mPhProj.DPOrient().SetDirIn(aModel.OriName());
    mPhProj.FinishInit();
    for (int aKIm=0 ; aKIm<2 ; aKIm++)
    {
        // Image names are stored as given to EpipRectification : current directory first, then project directory
        aStoredNames[aKIm] = aModel.ImName(aKIm+1);
        aNameIms[aKIm] = ExistFile(aStoredNames[aKIm]) ? aStoredNames[aKIm] : DirProject() + FileOfPath(aStoredNames[aKIm],false);
        MMVII_INTERNAL_ASSERT_User(ExistFile(aNameIms[aKIm]), eTyUEr::eOpenFile, "Image of the model not found : " + aStoredNames[aKIm]);
        aSIs[aKIm] = mPhProj.ReadSensor(FileOfPath(aStoredNames[aKIm],false /* Ok Not Exist*/),true/*DelAuto*/,false /* Not SVP*/);
    }

    const auto aCropOpts = EpipPairCropOpts(mCropMask,aModel.Map(1),aModel.Map(2),*aSIs[0],*aSIs[1],mOutDir,
                {EpipOutBaseName(mOutNamePat,aStoredNames[0],aStoredNames[1]),EpipOutBaseName(mOutNamePat,aStoredNames[1],aStoredNames[0])},true);

    for (int aKIm=0 ; aKIm<2 ; aKIm++)
    {
        StdOut() << Color::title << "* Image " << aKIm+1 << Color::end << std::endl;
        // RPC of the full frame saved by EpipRectification, next to the model : reused, else fitted again
        std::string aRPCFile;
        if (! mNoRPC)
        {
            const std::string aRPCName = aModel.RPCName(aKIm+1);
            if (! aRPCName.empty())
                aRPCFile = DirOfPath(mNameModel,false) + aRPCName;
            if (aRPCFile.empty() || ! ExistFile(aRPCFile))
            {
                MMVII_USER_WARNING("No saved RPC for image " + ToStr(aKIm+1) + " of the model, fitted again (" + aRPCName + ")");
                aRPCFile.clear();
            }
        }
        ResampleEpipImage(aCropOpts[aKIm], aNameIms[aKIm], aModel.Map(aKIm+1), mNoRPC ? nullptr : aSIs[aKIm], aRPCFile, *aInterp, false, mOutDir,
                          EpipOutBaseName(mOutNamePat,aStoredNames[aKIm],aStoredNames[1-aKIm]));
    }

    delete aInterp;
    return EXIT_SUCCESS;
}

int cAppli_EpipResampling::ExeTiles(const cEpipPairModel & aModel)
{
    MMVII_INTERNAL_ASSERT_User((mCropMask.mMaster==1) || (mCropMask.mMaster==2), eTyUEr::eUnClassedError, "Master must be 1 or 2");
    MMVII_INTERNAL_ASSERT_User((mSzTiles.x()>0) && (mSzTiles.y()>0), eTyUEr::eUnClassedError, "SzTiles must be positive");
    MMVII_INTERNAL_ASSERT_User((mSzOverL.x()>=0) && (mSzOverL.y()>=0) && (mSzOverL.x()<mSzTiles.x()) && (mSzOverL.y()<mSzTiles.y()),
        eTyUEr::eUnClassedError, "SzOverL must be non negative and smaller than SzTiles");

    // Region of the master epipolar frame : all of it, or the given crop
    cPt2di aR0(0,0), aR1 = aModel.Map(mCropMask.mMaster).EpipImSz();
    if (mCropMask.mHasCrop)
    {
        aR0 = mCropMask.mCropP0;
        aR1 = mCropMask.mCropP1;
    }
    const int aNbCol = (int)EpipTiles1D(aR0.x(),aR1.x(),mSzTiles.x(),mSzOverL.x()).size();
    const std::vector<cRect2> aTiles = EpipTiles(aR0,aR1,mSzTiles,mSzOverL);
    const int aNbRow = (int)aTiles.size() / aNbCol;
    StdOut() << "Tiling : " << aNbRow << " rows x " << aNbCol << " columns on [" << aR0 << " " << aR1 << "[" << std::endl;

    const std::string & aNameMaster = aModel.ImName(mCropMask.mMaster);
    const std::string & aNameSlave  = aModel.ImName(3-mCropMask.mMaster);

    // One call of this command per tile : the crop and the output name change, the other arguments are kept
    std::list<cParamCallSys> aListCom;
    std::vector<std::string> aInfoNames;
    for (size_t aKT=0 ; aKT<aTiles.size() ; aKT++)
    {
        const std::string aTilePat = EpipTileName(mOutNamePat,(int)aKT/aNbCol,(int)aKT%aNbCol,aNbRow,aNbCol);
        aInfoNames.push_back(EpipInfoName(mOutDir,EpipOutBaseName(aTilePat,aNameMaster,aNameSlave)));
        RemoveFile(aInfoNames.back(),SVP::Yes);   // a stale file must not hide a failure

        cParamCallSys aParam(cMMVII_Appli::FullBin(),Specs().Name(),mNameModel);
        for (size_t aKArg=3 ; aKArg<mArgv.size() ; aKArg++)   // [0] binary, [1] command, [2] the model (already given)
        {
            const std::string & aArg = mArgv[aKArg];
            bool aTilingArg = false;
            for (const std::string aName : {"SzTiles=","SzOverL=","CropP0=","CropP1=","OutName="})
                aTilingArg = aTilingArg || (aArg.rfind(aName,0)==0);
            if (! aTilingArg)
                aParam.AddArgs(aArg);
        }
        aParam.AddArgs("CropP0=" + cStrIO<cPt2di>::ToStr(aTiles[aKT].P0()));
        aParam.AddArgs("CropP1=" + cStrIO<cPt2di>::ToStr(aTiles[aKT].P1()));
        aParam.AddArgs("OutName=" + aTilePat);
        aParam.AddArgs(GIP_LevCall + "=" + ToStr(mLevelCall+1));
        aParam.AddArgs(GIP_KthCall + "=" + ToStr((int)aKT));
        aListCom.push_back(aParam);
    }
    const int aResult = ExeComParal(aListCom);

    // The index is built from the info files written by the tiles : only if all of them are there
    cEpipTilesInfo anIndex;
    anIndex.mNameModel = mNameModel;
    anIndex.mMaster = mCropMask.mMaster;
    anIndex.mSzTiles = mSzTiles;
    anIndex.mSzOverL = mSzOverL;
    anIndex.mRegion0 = aR0;
    anIndex.mRegion1 = aR1;
    cPt2dr aDisp(1e30,-1e30);
    std::string aFailed;
    for (size_t aKT=0 ; aKT<aTiles.size() ; aKT++)
    {
        if (! ExistFile(aInfoNames[aKT]))
        {
            aFailed += " (" + ToStr((int)aKT/aNbCol) + "," + ToStr((int)aKT%aNbCol) + ")";
            continue;
        }
        const cEpipCropInfo anInfo = cEpipCropInfo::FromFile(aInfoNames[aKT]);
        aDisp = cPt2dr(std::min(aDisp.x(),anInfo.mDispRange.x()),std::max(aDisp.y(),anInfo.mDispRange.y()));
        cEpipTileEntry anEntry;
        anEntry.mRow = (int)aKT/aNbCol;
        anEntry.mCol = (int)aKT%aNbCol;
        anEntry.mInfoFile = FileOfPath(aInfoNames[aKT],false);
        anIndex.mTiles.push_back(anEntry);
    }
    MMVII_INTERNAL_ASSERT_User(aFailed.empty(), eTyUEr::eUnClassedError, "Tiles (row,col) without result, no tile index written :" + aFailed);
    MMVII_INTERNAL_ASSERT_User(aResult==EXIT_SUCCESS, eTyUEr::eUnClassedError, "A command of the tiling failed, no tile index written");
    anIndex.mDispRange = aDisp;

    const std::string anIndexName = mOutDir + LastPrefix(EpipOutBaseName(mOutNamePat,aNameMaster,aNameSlave)) + ".EpipTiles." + GlobTaggedNameDefSerial();
    anIndex.ToFile(anIndexName);
    StdOut() << "Tile index : " << anIndexName << std::endl;
    return EXIT_SUCCESS;
}

cCollecSpecArg2007 & cAppli_EpipResampling::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return anArgObl
          << Arg2007(mNameModel,"epipolar model file of the image pair, as saved by EpipRectification's SaveModel=true",{eTA2007::FileTagged})
        ;
}

cCollecSpecArg2007 & cAppli_EpipResampling::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    anArgOpt
        << cHeaderSectionArg("Resampling")
        << AOpt2007(mInterpol,"Interpol","Interpolator", Append(cSpecOneArg2007::tAllSemPL{eTA2007::HDV},InterpolArgSem()))
        << AOpt2007(mCropMask.mCropP0,"CropP0","Crop : lower-left corner in the epipolar output coordinates (with CropP1 ; default : full frame)")
        << AOpt2007(mCropMask.mCropP1,"CropP1","Crop : upper-right corner in the epipolar output coordinates (with CropP0 ; default : full frame)")
        << AOpt2007(mCropMask.mMaster,"Master","Image (1 or 2) in which CropP0/CropP1 are given ; the crop of the other image is derived from the Z interval",{eTA2007::HDV})
        << cHeaderSectionArg("Tiling")
        << AOpt2007(mSzTiles,"SzTiles","Cut the pair in tiles of this size (master image pixels) and resample them in parallel (NbProc) ; region : master frame, or CropP0/CropP1",{eTA2007::HDV})
        << AOpt2007(mSzOverL,"SzOverL","Overlap between tiles (master image pixels)",{eTA2007::HDV})
        << cHeaderSectionArg("Output")
        << AOpt2007(mOutDir,"OutDir","Output directory (Default: VISU/" + Specs().Name()+")")
        << AOpt2007(mOutNamePat,"OutName","Output name pattern for images and other output files (%1 = image, %2 = other image ; must not be empty ; .tif is added if there is no extension)",{eTA2007::HDV})
        << AOpt2007(mNoRPC,"NoRPC","Don't write the RPC of the resampled images",{eTA2007::HDV})
        << AOpt2007(mCropMask.mMask,"Mask","Also write a 1-bit mask of the valid pixels of the resampled image (1 = maps inside the source image)",{eTA2007::HDV})
        << AOpt2007(mCropMask.mMaskName,"MaskName","Mask name pattern, $1 = output image name without extension ; implies Mask",{eTA2007::HDV})
        ;
    return anArgOpt;
}


std::vector<cOneHelpSampleCmp> cAppli_EpipResampling::Samples() const
{
    return
    {
        cOneHelpSampleCmp::Header("Resample the pair from the model saved by EpipRectification SaveModel=true"),
        {"MMVII EpipResampling VISU/EpipRectification/Epip_Im1_Im2.EpipModel.xml"},
        cOneHelpSampleCmp::Header("Crop (e.g. for dense matching) with masks"),
        {"MMVII EpipResampling VISU/EpipRectification/Epip_Im1_Im2.EpipModel.xml CropP0=[0,0] CropP1=[2000,1500] Mask=true"},
        cOneHelpSampleCmp::Header("Cut the pair in overlapping tiles resampled in parallel"),
        {"MMVII EpipResampling VISU/EpipRectification/Epip_Im1_Im2.EpipModel.xml SzTiles=[2000,2000] SzOverL=[100,100] NbProc=4"}
    };
}

tMMVII_UnikPApli Alloc_EpipResampling(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
    return tMMVII_UnikPApli(new cAppli_EpipResampling(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_EpipResampling
    (
        "EpipResampling",
        Alloc_EpipResampling,
        "Resample both images of an epipolar pair from a model saved by EpipRectification (SaveModel=true), with optional crop and validity mask",
        {eApF::ImProc},
        {eApDT::Orient,eApDT::Image},
        {eApDT::Orient,eApDT::Image},
        __FILE__
        );


}; // MMVII

