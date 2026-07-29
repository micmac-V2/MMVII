#include "BundleAdjustment.h"
#include "MMVII_Interpolators.h"
#include "MMVII_Tpl_Images.h"
#include "MMVII_ZBuffer.h"

namespace MMVII
{

//#define NUMPATCHDEBUG 3
//#define SCANSCANDEBUG 10
//#define EXPORTREPROJLIDAR

cBA_LidarBase::cBA_LidarBase(cPhotogrammetricProject * aPhProj,
                                     cMMVII_BundleAdj& aBA, double aSigma, const std::vector<std::string> & aInterp) :
    mPhProj     (aPhProj),
    mBA         (aBA),                                 // memorize the bundel adj class itself (access to optimizer)
    mParamInterpol (aInterp),
    mInterp     (cDiffInterpolator1D::AllocFromNames(mParamInterpol)),
    mEq         (nullptr),                             // equation of egalisation Lidar/Phgr
    mWFactor      (1/Square(aSigma)),
    mNbUsedPoints (0),
    mNbUsedObs (0),
    mNbUsableObs (0)
{

}

cBA_LidarBase::~cBA_LidarBase()
{
    delete mInterp;
}

//---------------------------------------------------

void cBA_LidarRaster::CreateZbuffers(cPhotogrammetricProject * aPhProj, const cMMVII_BundleAdj& aBA, bool aOnScans, bool aDebug)
{
    mMapZbuf.clear();
    // clear all patches unvisibility
    for (auto & aScanDataA: mVScans)
    {
        aScanDataA.mLidarRaster->TriangulateRegular(aPhProj->DirVisuAppli(),
                                                    1 + aScanDataA.mLidarRaster->InternalCalib()->F() / 100); // keep equivalent to a 0.01 rad step
        for (auto & aPatch:aScanDataA.mLPatches)
            aPatch.mHiddenOnImage.clear();
    }
    bool aZbufWithDist = true; // zbuffer is in dist or dz?
    std::vector<cSensorCamPC*> aVImages;

    // work on normal images (+ TSL only if aOnScans)
    for (const auto aPtrCam : aBA.VSCPC())
    {
        if (aOnScans || (dynamic_cast<cStaticLidar*>(aPtrCam)==nullptr))
        {
            aVImages.push_back(aPtrCam);
        }
    }

    std::string aVisuDir = aPhProj->DirVisuAppli()+aPhProj->DPOrient().DirIn() + "/";
    CreateDirectories(aVisuDir);
    for (auto & aCam: aVImages)
    {
        int aMarginInsideImage = 1;
        tREAL4 aDistTolerance = 0.2; // for overlapping walls with incorrect pose
        std::string aImName = aCam->NameImage();
        //std::cout<<"Visibility "<<aImName<<"\n";
        // create all z buffers
        for (auto & aScanDataA: mVScans)
        {
            auto &aScanA = aScanDataA.mLidarRaster;
            std::string aIndividualZbufFileName = "ZBuf-" + aScanA->NameImage()+"_on_"+aImName+".tif";

            if (ExistFile(aVisuDir + aIndividualZbufFileName))
            {
                cIm2D<tREAL4> aIndividualZbuf = cIm2D<tREAL4>::FromFile(aVisuDir + aIndividualZbufFileName);
                if (mMapZbuf.count(aImName)==0)
                    mMapZbuf.emplace(aImName, aIndividualZbuf);
                else {
                    // merge zbuffers
                    TransfoInPlace(mMapZbuf.at(aImName).DIm(), aIndividualZbuf.DIm(), [](auto &a, auto &b){ return std::max(a,b);} );
                }
                StdOut() << "Got zbuffer from: " << aVisuDir + aIndividualZbufFileName<<"\n";
                continue;
            }
            cMeshTri3DIterator aTriIt(aScanA->getTriangulation());
            if (aImName == aScanA->NameImage())
                continue;
            //ScopedTimer aTimer("Zbuffer");
            StdOut() << "Create zbuffer: " << aIndividualZbufFileName<<"\n";

            cSIMap_Ground2ImageAndProf aMapCamDepth(aCam);

            cSetVisibility aSetVis(aCam,aMarginInsideImage);

            double Infty =1e20;
            cPt2di aSzPix = aCam->SzPix();
            cBox3dr  aBox(cPt3dr(aMarginInsideImage,aMarginInsideImage,-Infty),
                         cPt3dr(aSzPix.x()-aMarginInsideImage,aSzPix.y()-aMarginInsideImage,Infty));
            cDataBoundedSet<tREAL8,3>  aSetCam(aBox);

            // full resolution zbuffer, to have the same coords as original image
            // TODO: smaller resolution
            cZBuffer aZBuf(aTriIt,aSetVis,aMapCamDepth,aSetCam,
                           1,aCam->PixelDomain().Sz());
            aZBuf.MakeZBuf(eZBufModeIter::ProjInit);
            aZBuf.ZBufIm().DIm().ToFile(aVisuDir + aIndividualZbufFileName);

            if (!aZBuf.IsOk())
            {
                StdOut() << "Warning! Zbuffer "<<aScanA->NameImage()<<" on image " << aImName <<" empty.\n";
                continue;
            }

            if (mMapZbuf.count(aImName)==0)
                mMapZbuf.emplace(aImName, aZBuf.ZBufIm());
            else {
                // merge zbuffers
                TransfoInPlace(mMapZbuf.at(aImName).DIm(), aZBuf.ZBufIm().DIm(), [](auto &a, auto &b){ return std::max(a,b);} );
            }
        }
        if (mMapZbuf.count(aImName)==0)
            return; // nothing more to do, we have no zbuffer info

        if (aDebug && (mMapZbuf.count(aImName)>0))
            mMapZbuf.at(aImName).DIm().ToFile(aPhProj->DirVisuAppli()
                                              +"ZBuf-total_image_"
                                              +aImName+".tif");

        // filter patches using all existing zbuffers
        for (auto & aScanDataA: mVScans)
        {
            //std::cout<<"test patchs of "<<aScanDataA.mLidarRaster->NameImage()<<"\n";
            for (auto & aPatch: aScanDataA.mLPatches)
            {
                const cPt2di & aPtScan = *aPatch.mLPatchesP.begin(); // center is 1st point
                cPt3dr aPtGround = aScanDataA.mLidarRaster->Image2Ground(aPtScan);
                cPt2dr aPtImage = aCam->Ground2Image(aPtGround);
                cPt3dr aPtCam3D = aCam->Pt_W2L(aPtGround);
                tREAL4 aDistWithTolerance = aZbufWithDist ?
                                            -(Norm2(aPtGround - aCam->Center()) - aDistTolerance) :
                                            -(aPtCam3D.z() - aDistTolerance);
                bool aIsUsablePt = true;

                //std::cout<<"  patch "<<aPatch.mId<<": "<<aPtScan<<" "<<aPtGround<<": in "<<aImName<<" "
                //          <<aPtImage<<": dist "<<Norm2(aPtGround - aCam->Center())
                //           <<" dz "<<aPtCam3D.z()<<"  final with tolerance: "<<aDistWithTolerance<<"\n";

                auto & aZbufIm = mMapZbuf.at(aImName).DIm();
                if (!aZbufIm.InsideBL(aPtImage))
                {
                    //std::cout<<"      out\n";
                    aIsUsablePt = false;

                } else {
                    auto aZbufVal = aZbufIm.GetVBL(aPtImage);
                    //std::cout<<"      zbuffer "<<aZbufVal<<"\n";
                    if ((aZbufVal > -1e9) && (aZbufVal > aDistWithTolerance))
                    {
                        //std::cout<<"      hidden\n";
                        aIsUsablePt = false;
                    }
                }
                //std::cout<<"        visible: "<<aIsUsablePt<<"\n";
                if (!aIsUsablePt)
                    aPatch.mHiddenOnImage.insert(aImName);
            }
        }
    }
}



//----------------------------------------------------------

cBA_LidarPhotogra::cBA_LidarPhotogra(cPhotogrammetricProject * aPhProj,
                                     cMMVII_BundleAdj& aBA, eImatchCrit aMode, double aSigma,
                                     const std::vector<std::string> & aInterp, bool aPertubate, int aNbPtsPerPatch) :
    cBA_LidarBase(aPhProj, aBA, aSigma, aInterp),
    mModeSim    (aMode),    // mode of matching
    mPertRad    (aPertubate),
    mNbPointByPatch (aNbPtsPerPatch)
{
    // read images before 1st iteration // TODO: read only images that may correspond to scans?
    StdOut() << "Read images...\n";
    for (const auto aPtrCam : aBA.VSCPC())
    {
        // do not read other TSLs
        if (dynamic_cast<cStaticLidar*>(aPtrCam))
            continue;
        auto & aImage = aPtrCam->LoadImage();
        if (mPertRad)
        {
            for (auto  aPix : aImage)
            {
                tREAL8 aMul =   (3+ sin(aPix.x()/70.0)) / 4.0;
                aImage.VI_SetV(aPix,aImage.VI_GetV(aPix)*aMul); // keep using ints?
            }
        }
    }

}

void cBA_LidarPhotogra::InitEq(bool aScanPoseUk)
{
    if (mModeSim==eImatchCrit::eDifRad)
        mEq = EqEqLidarImPonct (true,1,aScanPoseUk,true);
    else if (mModeSim==eImatchCrit::eCensus)
        mEq = EqEqLidarImCensus(true,1, aScanPoseUk,true);
    else if (mModeSim==eImatchCrit::eCorrel)
        mEq = EqEqLidarImCorrel(true,1, aScanPoseUk,true);
    else
    {
        MMVII_UnclasseUsEr("Bad enum for cBA_LidarPhotogra");
    }
}

cBA_LidarPhotograTri::cBA_LidarPhotograTri(cPhotogrammetricProject * aPhProj,
                                           cMMVII_BundleAdj& aBA,
                                           eImatchCrit aMode, const std::string & aPlyFile, double aSigma,
                                           const std::vector<std::string> & aInterp, bool aPertubate, int aNbPtsPerPatch) :
    cBA_LidarPhotogra(aPhProj, aBA, aMode, aSigma, aInterp, aPertubate, aNbPtsPerPatch), mTri(nullptr)
{
    InitEq(false);

    MMVII_INTERNAL_ASSERT_User((mModeSim!=eImatchCrit::eDifRad) || (mNbPointByPatch==1),
                               eTyUEr::eUnClassedError,"Only 1 point per patch in "+ToStr(eImatchCrit::eDifRad)+" mode");

    MMVII_INTERNAL_ASSERT_User(UCaseEqual(LastPostfix(aPlyFile),"ply"),
                               eTyUEr::eUnClassedError,"Lidar PLY file mandatory in triangulation mode, got \"" + aPlyFile + "\"");
    mTri = new cTriangulation3D<tREAL4>(aPlyFile);

    if (mModeSim!=eImatchCrit::eDifRad)
    {
        // cBox2dr aBox = BoxOfTri(mTri);
        cBox2dr aBox = mTri->Box2D();
        // estimate the distance for computing patching assuming a uniform  distributio,
        // Pi d^ 2  /NbByP = Surf / NbTot
        tREAL8 aDistMoy = std::sqrt(mNbPointByPatch *aBox.NbElem()/ (mTri->NbPts()*M_PI));
        tREAL8 aDistReject =  aDistMoy *1.5;
        mTri->MakePatches(mLPatchesI,aDistMoy,aDistReject,5);

        StdOut() << "Patches: DistReject=" << aDistReject
                 << " NbPts=" << mTri->NbPts()
                 << " NbPatch=" << mLPatchesI.size()
                 << "\n";
   }
}



cBA_LidarPhotograRaster::cBA_LidarPhotograRaster(cPhotogrammetricProject * aPhProj,
                                                 cMMVII_BundleAdj& aBA,
                                                 eImatchCrit aMode, const std::string & aPatScan, double aSigma,
                                                 const std::vector<std::string> & aInterp, double aScaleInit, double aScaleFinal,
                                                 double aThreshold, int aNbPtsPerPatch) :
    cBA_LidarPhotogra(aPhProj, aBA, aMode, aSigma, aInterp, false, aNbPtsPerPatch),
    mScaleInit(aScaleInit), mScaleFinal(aScaleFinal)
{
    InitEq(true);

    mThresholdInit = mThresholdFinal = aThreshold;

    tNameSelector aSel =   AllocRegex(aPatScan);
    for (const auto & aPtrCam : mBA.VSCPC())
    {
        cStaticLidar* aPtrScan = dynamic_cast<cStaticLidar*>(aPtrCam);
        if (aPtrScan && aSel.Match(aPtrScan->NameImage()) )
        {
            mBA.AddStaticLidar(aPtrScan);
            StdOut() << "Add Scan " << aPtrScan->NameImage() << "\n";
            mVScans.push_back({aPtrScan , {}});
            mIndexesScans[aPtrScan->NameImage()] = mVScans.size()-1;
        }
    }

    // Creation of the patches, choose a neigborhood around patch centers. TODO: adapt to images ground pixels size?
    if (mModeSim==eImatchCrit::eDifRad)
        mNbPointByPatch = 1;
    for (auto & aScanData: mVScans)
    {
        auto &aScan = aScanData.mLidarRaster;
        aScan->MakePatches(aScanData.mLPatches,aBA.VSCPC(),mNbPointByPatch,5,*mInterp);
        StdOut() << "Nb patches for " << aScan->NameImage() << ": " << aScanData.mLPatches.size() << "\n";
    }
}

cBA_LidarPhotogra::~cBA_LidarPhotogra()
{
}

cBA_LidarPhotograTri::~cBA_LidarPhotograTri()
{
    if (mTri) delete mTri;
}

cBA_LidarPhotograRaster::~cBA_LidarPhotograRaster()
{
}

void cBA_LidarPhotograTri::AddObs()
{
    mLastResidual.Reset();
    mNbUsedPoints = 0;
    mNbUsedObs = 0;
    mNbUsableObs = 0;
    cBasicWeighter<tREAL8> aWeighter(mWFactor);
    std::unordered_set<std::string> aNoHiddenPartComputed;
    if (mModeSim==eImatchCrit::eDifRad)
    {
        for (size_t aKP=0 ; aKP<mTri->NbPts() ; aKP++)
        {
            Add1Patch(aWeighter,{ToR(mTri->KthPts(aKP))},
                      "?",aNoHiddenPartComputed, aKP);
        }
    }
    else
    {
        int aNbPatch = 0;
        for (const auto& aPatchIndex : mLPatchesI)
        {
            std::vector<cPt3dr> aVP;
            for (const auto anInd : aPatchIndex)
                aVP.push_back(ToR(mTri->KthPts(anInd)));
            Add1Patch(aWeighter,aVP,"?",aNoHiddenPartComputed, aNbPatch);
            aNbPatch++;
        }
    }

    if (mLastResidual.SW() != 0)
       StdOut() << "  * Lid/Phr Residual Rad " << std::sqrt(mLastResidual.Average())
                 << " ("<<mNbUsedObs<<" obs, "<<mNbUsedPoints<<" points)\n";
    else
        StdOut() << "  * Lid/Phr: no obs\n";
}


//------------------------------------------------------

void cBA_LidarPhotograRaster::UpdateInterpolatorScale(const cMMVII_BundleAdj& aBA)
{
    tREAL4 aScale = aBA.NbMaxIter() < 2 ? mScaleFinal :
                     mScaleInit + (mScaleFinal - mScaleInit)*float(aBA.Iter())/(aBA.NbMaxIter()-1);
    std::cout << "up interpolator, scale="<<aScale<<"\n";

    if (mScaleInit == mScaleFinal)
        return;

    std::vector<std::string> aParamScaleInterpol = {"Scale", cStrIO<double>::ToStr(aScale), "1000"};
    aParamScaleInterpol.insert(aParamScaleInterpol.end(), mParamInterpol.begin(), mParamInterpol.end());
    delete mInterp;
    mInterp  = cDiffInterpolator1D::AllocFromNames(aParamScaleInterpol);
}


void cBA_LidarPhotograRaster::UpdateWeightersMap(const cMMVII_BundleAdj& aBA, double aWFactor)
{
    tREAL4 aTh = aBA.NbMaxIter() < 2 ? mThresholdFinal :
                     mThresholdInit + (mThresholdFinal - mThresholdInit)*float(aBA.Iter())/(aBA.NbMaxIter()-1);
    std::cout << "up weighters, th="<<aTh<<"\n";
    if (aTh>10000)
        aTh = -1;
    for (auto & aScanDataA: mVScans)
    {
        auto &aScanA = aScanDataA.mLidarRaster;
        tREAL8 aSigma = 1.; // TODO use image res for W? aScanA->Sigma() converted with incidence?
        mWeightersMap[aScanA->NameImage()].reset(new cStdWeighterResidual(sqrt(aWFactor)*aSigma, aTh / 2., aTh, 1));
    }
}


void cBA_LidarPhotograRaster::AddObs()
{
    mMapNbUsedPatches.clear();

    if (mBA.Iter()==0)
    {
        CreateZbuffers(mPhProj, mBA, false, true);
    }

    mLastResidual.Reset();
    mNbUsedPoints = 0;
    mNbUsedObs = 0;
    mNbUsableObs = 0;

    // update the interpolator and weighters map
    UpdateInterpolatorScale(mBA);
    UpdateWeightersMap(mBA, mWFactor);


    if (mModeSim==eImatchCrit::eDifRad)
    {
        for (auto & aScan : mVScans)
        {
            auto & aWeighter = mWeightersMap.at(aScan.mLidarRaster->NameImage());
            int aNbPatch = 0;
            for (const auto& aPatch : aScan.mLPatches)
            {
                Add1Patch(*aWeighter,
                          {aScan.mLidarRaster->Image2Ground(*aPatch.mLPatchesP.begin())},
                          aScan.mLidarRaster->NameImage(), aPatch.mHiddenOnImage, aNbPatch);
                aNbPatch++;
            }
        }
    }
    else
    {
        for (auto & aScan : mVScans)
        {
            auto & aWeighter = mWeightersMap.at(aScan.mLidarRaster->NameImage());
            int aNbPatch = 0;
            for (const auto& aPatch : aScan.mLPatches)
            {
                std::vector<cPt3dr> aVP;
                for (const auto aPt : aPatch.mLPatchesP)
                    aVP.push_back(aScan.mLidarRaster->Image2Ground(aPt));
                Add1Patch(*aWeighter,aVP,aScan.mLidarRaster->NameImage(), aPatch.mHiddenOnImage, aNbPatch);
                aNbPatch++;
            }
        }
    }

    if (mLastResidual.SW() != 0)
        StdOut() << "  * Lid/Phr Residual Rad " << std::sqrt(mLastResidual.Average())
                 << " ("<<mVScans.size()<<" scans, "<<mNbUsedObs<<" obs, "<<mNbUsedPoints<<" points)\n";
    else
        StdOut() << "  * Lid/Phr: no obs\n";

    if ((mBA.Iter()==0)||(mBA.Iter()==mBA.NbMaxIter()-1))
        for (const auto& [aCpl, aNb] : mMapNbUsedPatches)
            StdOut() <<  aCpl << ": " << aNb << " patches\n";
}


// ---------------------------------------------------------------------------

void cBA_LidarPhotograTri::SetVUkVObs
    (
     const cPt3dr&           aPGround,
     std::vector<int> *      aVIndUk,
     std::vector<tREAL8> &   aVObs,
     const cData1ImLidPhgr & aData,
     int                     aKPt
    )
{
    cSensorCamPC * aCam = mBA.VSCPC().at(aData.mKIm);  // extract the camera
    cPt3dr aPCam = aCam->Pt_W2L(aPGround);  // coordinate of point in image system
    tProjImAndGrad aPImGr = aCam->InternalCalib()->DiffGround2Im(aPCam); // compute proj & gradient

    // Vector of indexes of unknwons
    if (aVIndUk)
    {
       aCam->PushIndexes(*aVIndUk);       // add the unknowns [C,R] of the camera
    }


    // vector that will contains values of observation at this step
    aCam->Pose_WU().PushObs(aVObs,true);  // true because we transpose: we use W->C, which is the transposition of IJK : C->W

    aPGround.PushInStdVector(aVObs);   //
    aPCam.PushInStdVector(aVObs);
            
    aPImGr.mGradI.PushInStdVector(aVObs);  // Grad Proj/PCam
    aPImGr.mGradJ.PushInStdVector(aVObs);
            
    auto [aRad0,aGradIm] = aData.mVGr.at(aKPt);  // Radiom & grad
    aVObs.push_back(aRad0);
    aGradIm.PushInStdVector(aVObs);
}


void cBA_LidarPhotograRaster::SetVUkVObs
    (
        const cPt3dr&           aPGround,
        std::vector<int> *      aVIndUk,
        std::vector<tREAL8> &   aVObs,
        const cData1ImLidPhgr & aData,
        int                     aKPt
        )
{
    cStaticLidar * aScan = mBA.MapTSL().at(aData.mScanAName);
    cPt3dr aPScan = aScan->Pt_W2L(aPGround);  // coordinate of point in ground system
    cSensorCamPC * aCam = mBA.VSCPC().at(aData.mKIm);  // extract the camera
    cPt3dr aPCam = aCam->Pt_W2L(aPGround);  // coordinate of point in image system
    tProjImAndGrad aPImGr = aCam->InternalCalib()->DiffGround2Im(aPCam); // compute proj & gradient

    // Vector of indexes of unknwons
    if (aVIndUk)
    {
        aScan->PushIndexes(*aVIndUk);      // add the unknowns [C,R] of the scan
        aCam->PushIndexes(*aVIndUk);       // add the unknowns [C,R] of the camera
    }

    // vector that will contain values of observation at this step
    aScan->Pose_WU().PushObs(aVObs,false); // no transpose for scan
    aCam->Pose_WU().PushObs(aVObs,true);  // true because we transpose: we use W->C, which is the transposition of IJK : C->W

    aPScan.PushInStdVector(aVObs);   //
    aPCam.PushInStdVector(aVObs);

    aPImGr.mGradI.PushInStdVector(aVObs);  // Grad Proj/PCam
    aPImGr.mGradJ.PushInStdVector(aVObs);

    auto [aRad0,aGradIm] = aData.mVGr.at(aKPt);  // Radiom & grad
    aVObs.push_back(aRad0);
    aGradIm.PushInStdVector(aVObs);
}


std::pair<int,tREAL8> cBA_LidarPhotogra::AddPatchDifRad(const cBasicWeighter<tREAL8> &aWeighter,
     const std::vector<cPt3dr> & aVPatchPtGnd,
     const std::vector<cData1ImLidPhgr> &aVData,
     int aPatchNum)
{
     // read the solver now, because was not initialized at creation
     cResolSysNonLinear<tREAL8> *  aSys = mBA.Sys();

     cComputeStdDev<tREAL8>   aStdDev;    // compute the standard deviation of projected radiometry (indicator)

     cWeightAv<tREAL8,tREAL8> aWAv;       // compute average of image for radiom unknown
     for (const auto & aData : aVData)
     {
         tREAL8 aValIm = aData.mVGr.at(0).first;   // value of first/central pixel in this image
         aWAv.Add(1.0,aValIm);
         aStdDev.Add(1.0,aValIm);  // compute std deviation
     }


     cPt3dr    aPGround = aVPatchPtGnd.at(0);
     std::vector<tREAL8> aVTmpAvg{aWAv.Average()};  // vector for initializingz the temporay (here 1 = average)
     cSetIORSNL_SameTmp<tREAL8>  aStrSubst(aVTmpAvg); // structure for handling schurr eliminatio,
     // parse the data of the patch
     for (const auto & aData : aVData)
     {
         std::vector<int>       aVIndUk{-1}; // first one is a temporary (convention < 0)
         std::vector<tREAL8>    aVObs;
         SetVUkVObs (aPGround,&aVIndUk,aVObs,aData,0);
            
         // accumulate the equation involving the radiom
         aSys->R_AddEq2Subst(aStrSubst,mEq,aVIndUk,aVObs,aWeighter);
     }
     // do the substitution & add the equation reduced (Schurr complement)
     aSys->R_AddObsWithTmpUK(aStrSubst,mBA.CurLVMParam());
     return {aVData.size(), Square(aStdDev.StdDev(1e-5))};
}

std::pair<int, tREAL8> cBA_LidarPhotogra::AddPatchCensus(const cBasicWeighter<tREAL8> & aWeighter,
     const std::vector<cPt3dr> & aVPatchPtGnd,
     const std::vector<cData1ImLidPhgr> &aVData,
     int aPatchNum)
{
     // read the solver now, because was not initialized at creation
     cResolSysNonLinear<tREAL8> *  aSys = mBA.Sys();
     for (size_t aKPt=1; aKPt<aVPatchPtGnd.size() ; aKPt++)
     {
         // -------------- [1] Calculate the average ratio on all images --------------------
         cWeightAv<tREAL8,tREAL8> aAvRatio;  // stuct for averaging ratio
         for (const auto & aData : aVData)
         {
             tREAL8 aV0 = aData.mVGr.at(0).first;            // radiom of central pixel
             tREAL8 aVK = aData.mVGr.at(aKPt).first;         // radiom of neighbour
             aAvRatio.Add(1.0,NormalisedRatioPos(aV0,aVK)) ; // acumulate the ratio
         }
         std::vector<tREAL8> aVTmpAvg({aAvRatio.Average()});  // vector of value of temporary unknowns

         // -------------- [2] Add the observation --------------------
         cSetIORSNL_SameTmp<tREAL8>  aStrSubst(aVTmpAvg);  // structure for schur complement
         for (const auto & aData : aVData) // parse all the images
         {
             std::vector<int>  aVIndUk{-1} ;  // indexe of unknown
             std::vector<tREAL8>  aVObs;      // observation/context

             SetVUkVObs(aVPatchPtGnd.at(0)  ,&aVIndUk,aVObs,aData,0);            // add unkown AND observations
             SetVUkVObs(aVPatchPtGnd.at(aKPt),nullptr ,aVObs,aData,aKPt);        // add ONLY observations
             aSys->R_AddEq2Subst(aStrSubst,mEq,aVIndUk,aVObs,aWeighter); // add the equation in Schurr structure
         }
         // add all the equation to the system with Schurr's elimination
         aSys->R_AddObsWithTmpUK(aStrSubst,mBA.CurLVMParam());
     }
     return {aVData.size(), NAN};
}

std::pair<int, tREAL8> cBA_LidarPhotogra::AddPatchCorrel(const cBasicWeighter<tREAL8> &aWeighter,
     const std::vector<cPt3dr> & aVPatchPtGnd,
     const std::vector<cData1ImLidPhgr> &aVData,
     int aPatchNum)
{
     // read the solver now, because was not initialized at creation
     cResolSysNonLinear<tREAL8> *  aSys = mBA.Sys();
     // -------------- [1] Compute the normalized values --------------------
     size_t aNbPt = aVPatchPtGnd.size();
     //  vector that will store the normalized value (Avg=0, Sigma=1)
     cDenseVect<tREAL8>  aVMoy(aNbPt,eModeInitImage::eMIA_Null);

     //  memorize the radiometries of images as vector
     std::vector<cDenseVect<tREAL8>>  aListVRad;
     for (const auto & aData : aVData)
     {
         // change to vecor format
         cDenseVect<tREAL8> aV(aNbPt);
         for (size_t aK=0 ; aK< aNbPt ; aK++)
         {
             aV(aK)  = aData.mVGr.at(aK).first;
         }
         aListVRad.push_back(aV);
         cDenseVect<tREAL8> aV01 = NormalizeMoyNorm2(aV);  // noramlize value
         aVMoy += aV01;  //  accumulate in a vector
     }

     aVMoy *=  1/ tREAL8(aVData.size()); // make VMoy, average of normalized
     aVMoy =  NormalizeMoyNorm2(aVMoy);  // re normalized
       
     // -------------- [2] Intialize the temporary  --------------------

     /*  Say we have N points, M images,  tempory values will be stored "a la queue leu-leu" as :
               R1 .. RN  A0  B0 A1 B1 ... AM BM
             * where Ri are the unknown radiometry of the normalize patch
             * where Aj are the unkonw for tranfering radiom of image j to normalize patch such that

                 Ri =  Aj Imj(pij) + Bj

             Noting pij the projection of Pi in Imj
     */

     std::vector<tREAL8> aVTmp = aVMoy.ToStdVect(); // push first values of normalized patch
     int aK0Im = aVTmp.size();

     // push the initial values of Aj Bj
     for (const auto &  aVRad : aListVRad)
     {
         auto [A,B] =  LstSq_Fit_AxPBEqY(aVRad,aVMoy);  // solve  Ri = Aj Imj + Bj
         std::cout <<A<<" "<<B<<"\n";
         if (fabs(A)<1e-10)
             return {0,NAN}; // patch in a saturated area
         aVTmp.push_back(A); // add tmp unknown for Aj
         aVTmp.push_back(B); // add tmp unknown for Bj
     }

     cSetIORSNL_SameTmp<tREAL8>  aStrSubst(aVTmp); // structure for handling schurr eliminatio,

             // three structure for forcing conservation of normalizattion (Avg,Sigma) for VMoy
     std::vector<int> aVIndPt;       // indexe of unkown of norm radiom
     std::vector<tREAL8> aVFixAvg;   // vector for forcing average
     std::vector<tREAL8> aVFixVar;   // vector for forcing std dev

     // -------------- [3] Add the equation  --------------------


     for (int aKPt=0 ; aKPt <  (int) aNbPt ; aKPt++)  // parse all points
     {
         int aIndPt = -(1+aKPt);     // indexe of point are {-1,-2,....}
         aVIndPt.push_back(aIndPt);  // accumulat set of global indexe of unknown patch
         aVFixAvg.push_back(1.0);     //  Sum Rk = 0 => all weight = 1
         //  S(R+dR) ^ 2 =1   ;  S (2 R dR ) = 1 - S(R^2)  ; but S(R^2)=1 by construction ...
         aVFixVar.push_back(2*aVMoy(aKPt));

         for (int aKIm=0 ;  aKIm< (int) aVData.size() ; aKIm++)
         {
             int aIndIm = -(1+aK0Im+2*aKIm);  // compute indexe assumming "a la queue leu-leu"
             std::vector<int>       aVIndUk{aIndPt,aIndIm,aIndIm-1} ;  // indexes of 3 unknown
             std::vector<tREAL8>    aVObs;  // vector of observations
             SetVUkVObs (aVPatchPtGnd.at(aKPt),&aVIndUk,aVObs,aVData.at(aKIm),aKPt);  // read obs & global Uk
             aSys->R_AddEq2Subst(aStrSubst,mEq,aVIndUk,aVObs,aWeighter);  // add equation in tmp struct
         }
     }

     aStrSubst.AddOneLinearObs(aNbPt,aVIndPt,aVFixAvg,0.0);  // force average
     aStrSubst.AddOneLinearObs(aNbPt,aVIndPt,aVFixVar,0.0);  // force standard dev

     aSys->R_AddObsWithTmpUK(aStrSubst,mBA.CurLVMParam());
     return {aVData.size(),NAN};
}


void  cBA_LidarPhotogra::Add1Patch(const cBasicWeighter<tREAL8> &aWeighter,
                                  const std::vector<cPt3dr> & aVPatchPtGnd,
                                  const std::string & aScanName,
                                  const std::unordered_set<std::string> &aHiddenOnImage, int aPatchNum)
{
     std::vector<cData1ImLidPhgr> aVData; // for each image where patch is visible will store the data

     //  Parse all the image, we will select the images where all point of a patch are visible
     //std::cout<<"New patch\n";
     for (size_t aKIm=0 ; aKIm<mBA.VSCPC().size() ; aKIm++)
     {
          cSensorCamPC * aCam = mBA.VSCPC()[aKIm]; // extract cam
         if (dynamic_cast<cStaticLidar*>(aCam))
             continue;
          // 1st test: zbuffer visibility
          //std::cout<<"Im "<<aScanBData.mScanName<<" patch "<<aPatch.mId<<" vis "<<aPatch.mImVisible.at(aScanBData.mScanName)<<"\n";
          if (aHiddenOnImage.count(aCam->NameImage())>0)
              continue;

          auto & aGenDIm = aCam->LoadImage();
          if (aCam->IsVisible(aVPatchPtGnd.at(0))) // first test : is central point visible
          {
              cData1ImLidPhgr  aData; // data that will be filled
              aData.mScanAName = aScanName;
              aData.mKIm = aKIm;
              for (size_t aKPt=0 ; aKPt<aVPatchPtGnd.size() ; aKPt++) // parse the points of the patch
              {
                   cPt3dr aPGround = aVPatchPtGnd.at(aKPt);
                   if (aCam->IsVisible(aPGround))  // is the point visible in the camera
                   {
                        cPt2dr aPIm = mBA.VSCPC()[aKIm]->Ground2Image(aPGround); // extract the image  projection
                        if (aGenDIm.InsideInterpolator(*mInterp,aPIm,1.0))  // is it sufficiently inside
                        {
                            auto aVGr = aGenDIm.GetValueAndGradInterpol(*mInterp,aPIm); // extract pair Value/Grad of image
                            if ((aVGr.first==0)||(aVGr.first==255)) // refuse saturated pixels TODO improve criteria!
                                continue;
                            aData.mVGr.push_back(aVGr); // push it at end of stack
                        #ifdef NUMPATCHDEBUG
                            //if ((aPatchNum<=NUMPATCHDEBUG)&&(aData.mVGr.size()==1)) // show central point info
                            //    std::cout<<"Patch "<<aPatchNum<<" Gnd: "<<aPGround<< " on "<< mBA.VSCPC()[aKIm]->NameImage()<<": "<<aPIm<<" V= "<<aVGr.first<<"\n";
                        #endif
                        }
                   }
              }
              //  Does all the point of the patch were inside the image ?
              if (aData.mVGr.size() == aVPatchPtGnd.size())
              {
                  aVData.push_back(aData); // memorize the data for this image
              }
          }
     }
     //std::cout<<"\n";

     // if less than 2 images : nothing valuable to do
     if (aVData.size()<2) return;

#ifdef NUMPATCHDEBUG
     // debug patch
     int aPixSz = 15;
     int aSpaceSz = 1;
     if (aPatchNum<=NUMPATCHDEBUG)
     {
         for (const auto & aData : aVData)
         {
            if (aData.mVGr.size() != mNbPointByPatch)
                 continue; // not enougth data to export correctly
            int aW = sqrt(aData.mVGr.size());
            cRGBImage  aImDist8b(cPt2di(aW, aW)*(aPixSz+aSpaceSz)+cPt2di(aSpaceSz,aSpaceSz), cRGBImage::Gray128);
            int aI = 0;
            int aJ = 0;
            // make a vect of gray in correct order
            cDenseVect<tREAL8> aVGrOrdered(aData.mVGr.size());
            for (size_t i=0; i<aData.mVGr.size() ; ++i)
            {
                if (i==0)
                    aVGrOrdered(aData.mVGr.size()/2) = aData.mVGr[0].first; //center
                else if (i<=aData.mVGr.size()/2)
                    aVGrOrdered(i-1) = aData.mVGr[i].first;
                else
                    aVGrOrdered(i) = aData.mVGr[i].first;
            }

            cDenseVect<tREAL8> aVGrOrderedNormed = NormalizeMoyVar(aVGrOrdered);  // noramlize value
            auto aMin = aVGrOrderedNormed.Min();
            auto aMax = aVGrOrderedNormed.Max();

            for (size_t i=0; i<aData.mVGr.size(); ++i)
            {
                auto aV = (aVGrOrderedNormed(i)-aMin)/(aMax-aMin)*255;
                aImDist8b.FillRectangle(cPt3di(aV,aV,aV),
                                        cPt2di(aI, aJ)*(aPixSz+aSpaceSz)+cPt2di(aSpaceSz,aSpaceSz),
                                        cPt2di(aI+1, aJ+1)*(aPixSz+aSpaceSz),
                                        cPt3dr(0.,0.,0.));
                //aImDist8b.RawSetPoint(cPt2di(aI, aJ)*(aPixSz+aSpaceSz), aV, aV, aV);
                aI++;
                if (aI==aW)
                {
                    aJ++;
                    aI = 0;
                }
            }
            std::string aPath = mPhProj->DirVisuAppli() + "patch_" + ToStr(aPatchNum,3) + "_iter" + ToStr(mBA.Iter(),1)
                                + "_" + mBA.VSCPC()[aData.mKIm]->NameImage() + ".png";
            aImDist8b.ToFile(aPath);
         }
     }
#endif

     mNbUsedPoints++;

     // accumlulate for computing average of deviation
     // mLastResidual.Add(1.0,  (aStdDev.StdDev(1e-5) *aVData.size()) / (aVData.size()-1.0));
     // mLastResidual.Add(1.0,  (aStdDev.StdDev(1e-5) ) );
     //mLastResidual.Add(aVData.size(),  Square(aStdDev.StdDev(1e-5) ) );

     int aNb = 0;
     tREAL8 aRes2 = 0.;
     if (mModeSim==eImatchCrit::eDifRad)
     {
         std::tie(aNb, aRes2) = AddPatchDifRad(aWeighter,aVPatchPtGnd,aVData, aPatchNum);
     }
     else if (mModeSim==eImatchCrit::eCensus)
     {
         std::tie(aNb, aRes2) = AddPatchCensus(aWeighter,aVPatchPtGnd,aVData, aPatchNum);
     }
     else if (mModeSim==eImatchCrit::eCorrel)
     {
         std::tie(aNb, aRes2) = AddPatchCorrel(aWeighter,aVPatchPtGnd,aVData, aPatchNum);
     }
     mNbUsedObs+=aNb;
     mLastResidual.Add(aNb,  aRes2 );
}

//-------------------------------------------------------------


std::pair<int, tREAL8> cBA_LidarPhotograRaster::AddPatchCorrel(const cBasicWeighter<tREAL8> &aWeighter,
     const std::vector<cPt3dr> & aVPatchPtGnd,
     const std::vector<cData1ImLidPhgr> &aVData,
     int aPatchNum)
{
    // read the solver now, because was not initialized at creation
    cResolSysNonLinear<tREAL8> *  aSys = mBA.Sys();
    // -------------- [1] Compute the normalized values --------------------
    size_t aNbPt = aVPatchPtGnd.size();
    //  vector that will store the normalized value (Avg=0, Sigma=1)
    cDenseVect<tREAL8>  aVMedian(aNbPt,eModeInitImage::eMIA_Null);

    //  memorize the radiometries of images as vector
    std::vector<cDenseVect<tREAL8>>  aListVRad;
    // memorize normalized radiometries for median
    std::vector<cDenseVect<tREAL8>>  aListVRadNorm;
    for (const auto & aData : aVData)
    {
        // change to vecor format
        cDenseVect<tREAL8> aV(aNbPt);
        for (size_t aK=0 ; aK< aNbPt ; aK++)
        {
            aV(aK)  = aData.mVGr.at(aK).first;
        }
        aListVRad.push_back(aV);
        aListVRadNorm.push_back(NormalizeMoyNorm2(aV));  // normalize value

    }
    for (size_t aK=0 ; aK< aNbPt ; aK++)
    {
        std::vector<tREAL8> aVV(aVData.size());
        for (size_t aJ=0 ; aJ< aVData.size() ; aJ++)
        {
            aVV[aJ] = aListVRadNorm[aJ](aK);
        }
        aVMedian(aK) = NonConstMediane(aVV);
    }
    aVMedian =  NormalizeMoyNorm2(aVMedian);  // re normalized

    // -------------- [2] Intialize the temporary  --------------------

    /*  Say we have N points, M images,  tempory values will be stored "a la queue leu-leu" as :
               R1 .. RN  A0  B0 A1 B1 ... AM BM
             * where Ri are the unknown radiometry of the normalize patch
             * where Aj are the unkonw for tranfering radiom of image j to normalize patch such that

                 Ri =  Aj Imj(pij) + Bj

             Noting pij the projection of Pi in Imj
     */

    std::vector<tREAL8> aVTmp = aVMedian.ToStdVect(); // push first values of normalized patch
    int aK0Im = aVTmp.size();

    // push the initial values of Aj Bj
    std::vector<int> aVIndexUsedImages;
    int aNumIm = -1;
    tREAL8 aMeanRes2 = 0.;
    for (const auto &  aVRad : aListVRad)
    {
        aNumIm++;
        auto [A,B] =  LstSq_Fit_AxPBEqY(aVRad,aVMedian);  // solve  Ri = Aj Imj + Bj
        // get residuals
        cDenseVect<tREAL8> aVect1(aNbPt,eModeInitImage::eMIA_V1);
        auto aRes = sqrt((A * aVRad + aVect1*B - aVMedian).SqL2Norm(false)); //quadratic mean residual for stddev=1 images
        // zncc = 1-aRes/2 ?
        auto aW = aWeighter.WeightOfResidual({aRes})[0];
#ifdef NUMPATCHDEBUG
        if (aPatchNum<=NUMPATCHDEBUG)
            StdOut() << "patch " << aPatchNum << " im " << aNumIm << " A=" << A << " B=" << B << " res=" << aRes << " W=" << aW << "\n";
#endif
        if (aW==0)
        {
#ifdef NUMPATCHDEBUG
            if (aPatchNum<=NUMPATCHDEBUG)
                StdOut() << "patch image rejected\n";
#endif
            continue; // patch too far
        }
        if (fabs(A)<1e-10)
            return {0,0.}; // patch in a saturated area
        aVTmp.push_back(A); // add tmp unknown for Aj
        aVTmp.push_back(B); // add tmp unknown for Bj
        aVIndexUsedImages.push_back(aNumIm);
        aMeanRes2 += aRes*aRes;
    }
    if (aVIndexUsedImages.size()<2)
        return {0,0.}; // this patch does not have enought suitable images

    for (auto & aData: aVData)
    {
        auto aCplId = std::pair(aData.mScanAName, mBA.VSCPC().at(aData.mKIm)->NameImage());
        if (mMapNbUsedPatches.count(aCplId)==0)
            mMapNbUsedPatches[aCplId] = 1;
        else
            mMapNbUsedPatches[aCplId]++;
    }

    // use the same weight for each point eq
    auto aRes = sqrt(aMeanRes2/aVIndexUsedImages.size());
    auto aW = aWeighter.WeightOfResidual({aRes})[0];
#ifdef NUMPATCHDEBUG
    if (aPatchNum<=NUMPATCHDEBUG)
        StdOut() << "patch " << aPatchNum << " res=" << aRes << " W=" << aW << "\n";
#endif
    cSetIORSNL_SameTmp<tREAL8>  aStrSubst(aVTmp); // structure for handling schurr eliminatio,

        // three structure for forcing conservation of normalizattion (Avg,Sigma) for VMoy
        std::vector<int> aVIndPt;       // indexe of unkown of norm radiom
    std::vector<tREAL8> aVFixAvg;   // vector for forcing average
    std::vector<tREAL8> aVFixVar;   // vector for forcing std dev

    // -------------- [3] Add the equation  --------------------


    for (int aKPt=0 ; aKPt <  (int) aNbPt ; aKPt++)  // parse all points
    {
        int aIndPt = -(1+aKPt);     // indexe of point are {-1,-2,....}
        aVIndPt.push_back(aIndPt);  // accumulat set of global indexe of unknown patch
        aVFixAvg.push_back(1.0);     //  Sum Rk = 0 => all weight = 1
        //  S(R+dR) ^ 2 = 1  ;  S (2 R dR ) = 1 - S(R^2)  ; but S(R^2)=1 by construction ...
        aVFixVar.push_back(2*aVMedian(aKPt));

        int aNumUsedIm = 0; // for uk index
        for (auto aKIm: aVIndexUsedImages)
        {
            int aIndIm = -(1+aK0Im+2*aNumUsedIm);  // compute indexe assumming "a la queue leu-leu"
            std::vector<int>       aVIndUk{aIndPt,aIndIm,aIndIm-1} ;  // indexes of 3 unknown
            std::vector<tREAL8>    aVObs;  // vector of observations
            SetVUkVObs (aVPatchPtGnd.at(aKPt),&aVIndUk,aVObs,aVData.at(aKIm),aKPt);  // read obs & global Uk
            //std::cout<<"ind: ";
            //for (auto &v: aVIndUk)
            //    std::cout<<v<<" ";
            //std::cout<<"\n";
            aSys->R_AddEq2Subst(aStrSubst,mEq,aVIndUk,aVObs,aW);  // add equation in tmp struct
            aNumUsedIm++;
        }
    }

    // TODO: check weights with patch size
    auto aW0 = aWeighter.WeightOfResidual({0.})[0];

    aStrSubst.AddOneLinearObs(aW0*aNbPt,aVIndPt,aVFixAvg,0.0);  // force average
    aStrSubst.AddOneLinearObs(aW0*aNbPt,aVIndPt,aVFixVar,0.0);  // force standard dev

    aMeanRes2 = aMeanRes2/(aVIndexUsedImages.size()*aVIndexUsedImages.size());
    aSys->R_AddObsWithTmpUK(aStrSubst,mBA.CurLVMParam());
    return {aVIndexUsedImages.size(),aMeanRes2};
}


//-------------------------------------------------------------

cBA_LidarLidarRaster::cBA_LidarLidarRaster(cPhotogrammetricProject * aPhProj,
                                           cMMVII_BundleAdj& aBA, const std::string & aPatScan, double aSigma,
                                           double aThresholdInit, double aThresholdFinal,
                                           double aNormalTolDeg, const std::vector<std::string> & aInterp) :
    cBA_LidarBase(aPhProj, aBA, aSigma, aInterp)
{
    mEq = EqEqLidarLidar (true,1,true);

    mThresholdInit = (aThresholdInit<0) ? INFINITY : aThresholdInit;
    mThresholdFinal = (aThresholdFinal<0) ? INFINITY : aThresholdFinal;

    MMVII_INTERNAL_ASSERT_User((aNormalTolDeg>=0) && (aNormalTolDeg<=180),
                               eTyUEr::eBadOptParam,"Normal tolerance must be inside [0,180], got "+ToStr(aNormalTolDeg));
    mNormalDiffMinCos = cos(aNormalTolDeg*M_PI/180.);

    //read scans files from directory corresponding to pattern in aPatScan
    tNameSelector aSel =   AllocRegex(aPatScan);
    for (const auto & aPtrCam : mBA.VSCPC())
    {
        cStaticLidar* aPtrScan = dynamic_cast<cStaticLidar*>(aPtrCam);
        if (aPtrScan && aSel.Match(aPtrScan->NameImage()) )
        {
            mBA.AddStaticLidar(aPtrScan);
            StdOut() << "Add Scan " << aPtrScan->NameImage() << "\n";
            mVScans.push_back({aPtrScan , {}});
            mIndexesScans[aPtrScan->NameImage()] = mVScans.size()-1;
        }
    }

    MMVII_INTERNAL_ASSERT_User(!mVScans.empty(),
                               eTyUEr::eBadFileSetName,"No TSL found!");

    // Creation of the patches, here just center point
    for (auto & aScanData: mVScans)
    {
        auto & aScan = aScanData.mLidarRaster;
        aScan->MakePatches(aScanData.mLPatches,aBA.VSCPC(),1,5, *mInterp);
        StdOut() << "Nb patches for " << aScan->NameImage()<< ": " << aScanData.mLPatches.size() << "\n";

        //for (auto &aTestRasterPoint: {cPt2di(10672,2238), cPt2di(2552,2121) })
        //    StdOut() << "Test " << aScanData.mScanName << " " << aTestRasterPoint << ": "
        //             << aScanData.mLidarRaster->Image2Ground(aTestRasterPoint) <<"\n";
    }
}

cBA_LidarLidarRaster::~cBA_LidarLidarRaster()
{
}


void cBA_LidarLidarRaster::UpdateWeightersMap(const cMMVII_BundleAdj& aBA, double aWFactor)
{
    tREAL4 aTh = aBA.NbMaxIter() < 2 ? mThresholdFinal :
                     mThresholdInit + (mThresholdFinal - mThresholdInit)*float(aBA.Iter())/(aBA.NbMaxIter()-1);
    //std::cout << "up weighters, th="<<aTh<<"\n";
    if (aTh>10000)
        aTh = -1;
    for (auto & aScanDataA: mVScans)
    {
        auto &aScanA = aScanDataA.mLidarRaster;
        for (auto & aScanDataB: mVScans)
        {
            auto &aScanB = aScanDataB.mLidarRaster;
            tREAL8 aSigmaAB = sqrt(aScanA->Sigma()*aScanA->Sigma()
                                   +aScanB->Sigma()*aScanB->Sigma());
            mWeightersMap[aScanA->NameImage()+"-"+aScanB->NameImage()].reset(
                //= cStdWeighterResidual(sqrt(aWFactor)*aSigmaAB, aTh / 30., aTh, 1); // simulate least L1 with 1/31 of weight before exclusion
                //= cStdWeighterResidual(sqrt(aWFactor)*aSigmaAB, aTh / 9., aTh, 1); // simulate least L1 with 1/10 of weight before exclusion
                //new cStdWeighterResidual(sqrt(aWFactor)*aSigmaAB, -1, aTh, 1) // least squares better for final convergence since exclusion methods are efficient
                new cLinearWeighterResidual(sqrt(aWFactor)*aSigmaAB, aTh, aTh*10)
                );
        }
    }
}


void cBA_LidarLidarRaster::AddObs()
{
    mMapNbUsedPatches.clear();

    if (mBA.Iter()>=0)
    {
        //CreateZbuffers(mPhProj, mBA, true, true); // useless for lidarlidar
    }

#ifdef SCANSCANSHOWPATCHES
    if (mBA.Iter()==mBA.NbMaxIter()-1)
    {
        for (const auto& aScanA : mVScans)
            for (const auto& aScanB : mVScans)
                if (aScanA.mLidarRaster!=aScanB.mLidarRaster)
                    mMapPatchesRasters.try_emplace(aScanA.mLidarRaster->NameImage()+"_to_"+aScanB.mLidarRaster->NameImage(),
                                                   aScanA.mLidarRaster->InternalCalib()->SzPix()/SCANSCANSHOWPATCHES + cPt2di(1,1),
                                                   nullptr,eModeInitImage::eMIA_Null);
    }
#endif

    mLastResidual.Reset();
    mNbUsedPoints = 0;
    mNbUsedObs = 0;
    mNbUsableObs = 0;

    // update the weighters map
    UpdateWeightersMap(mBA, mWFactor);

    for (auto & aScan : mVScans)
    {
#ifdef SCANSCANDEBUG
        cIm2D<tREAL4> aResImage(aScan.mLidarRaster->InternalCalib()->SzPix(),0,eModeInitImage::eMIA_Null);
        auto & aResImageData = aResImage.DIm();
        int aPtSize = 1 + aScan.mLidarRaster->PixelDomain().Sz().x()/1000;
        if (mBA.Iter()%SCANSCANDEBUG==0)
        {
            for (int y=0; y<aScan.mLidarRaster->InternalCalib()->SzPix().y();++y)
                for (int x=0; x<aScan.mLidarRaster->InternalCalib()->SzPix().x();++x)
                    aResImageData.SetV({x,y}, 999);
        }
#endif

        for (const auto& aPatch : aScan.mLPatches)
        {
            //if (*aPatch.begin()==cPt2di(4278, 2245)) //10677, 2481
            //    std::cout<<"!\n";
            [[maybe_unused]] auto aMinRes = Add1Patch(aPatch, aScan);
#ifdef SCANSCANDEBUG
            if (mBA.Iter()%SCANSCANDEBUG==0)
            {
                auto aC = *aPatch.mLPatchesP.begin();
                for (int y=aC.y()-aPtSize+1; y<aC.y()+aPtSize;++y)
                    for (int x=aC.x()-aPtSize+1; x<aC.x()+aPtSize;++x)
                        aResImageData.SetVTruncIfInside({x,y}, std::isnan(aMinRes)?999:aMinRes);
            }
#endif
        }

#ifdef SCANSCANDEBUG
        if (mBA.Iter()%SCANSCANDEBUG==0)
        {
            std::string aPath = mPhProj->DirVisuAppli() + aScan.mLidarRaster->NameImage() + "_iter_" + ToStr(mBA.Iter())+ ".tif";
            aResImageData.ToFile(aPath, {"COMPRESS=DEFLATE"});
        }
#endif

    }
    if (mLastResidual.SW() != 0)
    {
        StdOut() << "  * Lid/Lid Residual dist " << std::sqrt(mLastResidual.Average())
                 << "m ("<<mVScans.size()<<" scans, "<<mNbUsedObs<<" obs="<< 100.*mNbUsedObs/mNbUsableObs
                 <<"%, "<<mNbUsedPoints<<" points)\n";
        //for (auto & aScan:mVScans)
        //    StdOut() << aScan.mLidarRaster->NameImage()<< " " << aScan.mLidarRaster->Center().x() <<
        //         std::setprecision(10) << " " << aScan.mLidarRaster->Center().y()<< " " << aScan.mLidarRaster->Center().z() << "\n";
    }
    else
        StdOut() << "  * Lid/Lid: no obs\n";

    if ((mBA.Iter()==0)||(mBA.Iter()==mBA.NbMaxIter()-1))
    {
        //for (const auto& [aCpl, aNb] : mMapNbUsedPatches)
        //    StdOut() <<  aCpl << ": " << aNb << " patches\n";

        tDMatR aNbPatchesMat(mVScans.size(),mVScans.size(),eModeInitImage::eMIA_Null);
        for (const auto& [aCpl, aNb] : mMapNbUsedPatches)
            aNbPatchesMat.SetElem(mIndexesScans.at(aCpl.first),mIndexesScans.at(aCpl.second),aNb);
        //StdOut() << aNbPatchesMat <<"\n";

        std::vector<std::string> aVNames;
        for (const auto& aScan : mVScans)
            aVNames.push_back(SplitString(aScan.mLidarRaster->NameImage(),".").at(0));
        StdOut() << "Patches visibility (col to row):\n";

        ShowMatrixWithNames(aNbPatchesMat, aVNames, aVNames);
    }

#ifdef SCANSCANSHOWPATCHES
    if (mBA.Iter()==mBA.NbMaxIter()-1)
    {
        for (const auto& [aCpl, aIm] : mMapPatchesRasters)
        {
            std::string aPath = mPhProj->DirVisuAppli() + "Patches_iter_" + ToStr(mBA.Iter()) + "_" + aCpl + ".tif";
            cIm2D<tREAL4> aZoomed = aIm.EnlargeInt(SCANSCANSHOWPATCHES);
            aZoomed.DIm().ToFile(aPath, {"COMPRESS=DEFLATE"});
        }
    }
#endif

#ifdef EXPORTREPROJLIDAR
    if (mBA.Iter()==mBA.NbMaxIter()-1)
    {
        // reproject intensity for scans with common patches
        for (const auto& [aCpl, aNb] : mMapNbUsedPatches)
        {
            std::string aPath = mPhProj->DirVisuAppli() + "Reproj_on_" + aCpl.first + "_intensity_from_" + aCpl.second + ".tif";
            cStaticLidar* aScanA =  mBA.MapTSL().at(aCpl.first);
            cStaticLidar* aScanB =  mBA.MapTSL().at(aCpl.second);
            auto aProjection = aScanA->projectIntensityFrom(*aScanB);
            aProjection.DIm().ToFile(aPath, {"COMPRESS=DEFLATE"});
        }
    }
#endif
}


void cBA_LidarLidarRaster::SetVUkVObs
    (const cPt3dr&           aPGround,
     std::vector<int> *      aVIndUk,
     std::vector<tREAL8> &   aVObs,
     const cData1ImLidPhgr & aData,
     int                     aKPt
     )
{
    cStaticLidar * aScanA = mBA.MapTSL().at(aData.mScanAName);
    cPt3dr aPScanA = aScanA->Pt_W2L(aPGround);  // coordinate of point in ground system
    cStaticLidar * aScanB = mBA.MapTSL().at(aData.mScanBName);
    cPt3dr aPScanB0 = aScanB->Pt_W2L(aPGround);  // coordinate of point in image system
    tProjImAndGrad aPImGr = aScanB->InternalCalib()->DiffGround2Im(aPScanB0); // compute proj & gradient

    // Vector of indexes of unknwons
    if (aVIndUk)
    {
        aScanA->PushIndexes(*aVIndUk);      // add the unknowns [C,R] of the scan
        aScanB->PushIndexes(*aVIndUk);       // add the unknowns [C,R] of the camera
    }

    // vector that will contain values of observation at this step
    aScanA->Pose_WU().PushObs(aVObs,false); // no transpose for scan
    aScanB->Pose_WU().PushObs(aVObs,true);  // true because we transpose: we use W->C, which is the transposition of IJK : C->W

    aPScanA.PushInStdVector(aVObs);   //
    aPScanB0.PushInStdVector(aVObs);

    aPImGr.mGradI.PushInStdVector(aVObs);  // Grad Proj/PCam
    aPImGr.mGradJ.PushInStdVector(aVObs);

    auto [aRad0,aGradIm] = aData.mVGr.at(aKPt);  // Radiom & grad
    aVObs.push_back(aRad0);
    aGradIm.PushInStdVector(aVObs);
}

tREAL8 cBA_LidarLidarRaster::Add1Patch(const cLidarRasterPatch &aPatch, const cStaticLidarBAData & aScanAData)
{
    auto & aScanA = aScanAData.mLidarRaster;
    cPt3dr aPGround = aScanA->Image2Ground(*aPatch.mLPatchesP.begin());
    std::vector<cData1ImLidPhgr> aVData; // for each image where patch is visible will store the data
    cWeightAv<tREAL8>   aAvgRes;    // compute average residual
    tREAL8 aMinResidual = INFINITY;
    cPt3dr aNormalGndA = aScanA->Pose().Rot().Value(aPatch.mNormalInstr);

    #ifdef SCANSCANDEBUG
    std::cout<<"ScanA: "<<aScanA->NameImage()<<" Patch "<<aPatch.mId<<": "<<*aPatch.mLPatchesP.begin()<<" -> Gnd: "<<aPGround<<"\n";
    #endif
    //  Parse all the scans B, we will select the ones where the patch is visible
    for (auto & aScanBData: mVScans)
    {
        auto & aScanB = aScanBData.mLidarRaster;
        auto & aWeighter = mWeightersMap.at(aScanA->NameImage()+"-"+aScanB->NameImage());
        if (aScanB->NameImage()==aScanA->NameImage())
            continue; // no obs on the same scan

        // 1st test: zbuffer visibility
        #ifdef SCANSCANDEBUG
        std::cout<<"On scan B "<<aScanB->NameImage()<<": ";
        #endif
        if (aPatch.mHiddenOnImage.count(aScanB->NameImage())>0)
        {
            #ifdef SCANSCANDEBUG
            std::cout<<" hidden\n";
            #endif

            #ifdef SCANSCANSHOWPATCHES
            if (mBA.Iter()==mBA.NbMaxIter()-1)
                mMapPatchesRasters.at(aScanA->NameImage()+"_to_"+aScanB->NameImage()).DIm().SetV(
                    aPatch.mLPatchesP[0]/SCANSCANSHOWPATCHES, -30); // hidden zbuffer
            #endif
            continue;
        }
        cDataGenUnTypedIm<2> & aGenDImDist = aScanB->getRasterDistance();

        // recheck if central point visible, TODO: remove, aPatch.mImVisible should be enought
        if (aScanB->IsVisible(aPGround))
        {
            cData1ImLidPhgr  aData; // data that will be filled
            aData.mScanAName = aScanA->NameImage();
            aData.mScanBName = aScanB->NameImage();
            cPt2dr aPIm = aScanB->Ground2Image(aPGround); // extract the image  projection
            #ifdef SCANSCANDEBUG
            std::cout<<" projection :"<<aPIm<<"\n";
            #endif
            tREAL8 aDist = Norm2(aPGround-aScanB->Center());
            if (aScanB->IsMaskedPoint(aPIm))
            {
                #ifdef SCANSCANDEBUG
                std::cout<<" masked point\n";
                #endif
                #ifdef SCANSCANSHOWPATCHES
                if (mBA.Iter()==mBA.NbMaxIter()-1)
                    mMapPatchesRasters.at(aScanA->NameImage()+"_to_"+aScanB->NameImage()).DIm().SetV(
                        aPatch.mLPatchesP[0]/SCANSCANSHOWPATCHES, -10); // masked
                #endif
                continue;
            }
            if (aGenDImDist.InsideInterpolator(*mInterp,aPIm,1.0))  // is it sufficiently inside
            {
                auto aVGr = aGenDImDist.GetValueAndGradInterpol(*aScanB->getLineraInterpolator(),aPIm); // extract pair Value/Grad of image

                aData.mVGr = {aVGr};
                #ifdef SCANSCANDEBUG
                std::cout<<aPIm<<" GetValueAndGradInterpol: "<<aVGr.first<<" "<<aVGr.second.x()*1940<<" "<<aVGr.second.y()*1940<<"\n";
                #endif

                tREAL8 aValIm = aData.mVGr.at(0).first;   // value of first/central pixel in this image
                tREAL8 aResidual = aValIm-aDist;

                if (fabs(aResidual)<std::max(0.1,mThresholdInit*20)) // suppose that 10cm is always an error
                    mNbUsableObs++;

                if (fabs(aResidual)<fabs(aMinResidual))
                    aMinResidual = aResidual;
                //StdOut() << "UUU  " << ((cStdWeighterResidual*)&aWeighter)->SingleWOfResidual(std::vector<tREAL8>{aResidual})
                //         << " " << aWeighter->WeightOfResidual({aResidual})[0] << std::endl;
                if (aWeighter->WeightOfResidual({aResidual})[0]==0.0)
                {
                    #ifdef SCANSCANDEBUG
                    std::cout<<"removed W\n";
                    #endif
                    #ifdef SCANSCANSHOWPATCHES
                    if (mBA.Iter()==mBA.NbMaxIter()-1)
                        mMapPatchesRasters.at(aScanA->NameImage()+"_to_"+aScanB->NameImage()).DIm().SetV(
                            aPatch.mLPatchesP[0]/SCANSCANSHOWPATCHES, 10 + fabs(aResidual)); // rejected for residual
                    #endif
                    continue;
                }
                cPt3dr aNormalInstrB = aScanB->Image2NormalInstr(aPIm, *mInterp);
                cPt3dr aNormalGndB = aScanB->Pose().Rot().Value(aNormalInstrB);
                if (Scal(aNormalGndA,aNormalGndB)<mNormalDiffMinCos)
                {
                    #ifdef SCANSCANDEBUG
                    std::cout<<"Removed "<<aPatch.mLPatchesP[0]<<" due to normals: "<<aNormalGndA<<" "<<aNormalGndB
                             <<" "<<acos(Scal(aNormalGndA,aNormalGndB))*180/M_PI<<"deg\n";
                    #endif
                    #ifdef SCANSCANSHOWPATCHES
                    if (mBA.Iter()==mBA.NbMaxIter()-1)
                        mMapPatchesRasters.at(aScanA->NameImage()+"_to_"+aScanB->NameImage()).DIm().SetV(
                            aPatch.mLPatchesP[0]/SCANSCANSHOWPATCHES, -1000 - acos(Scal(aNormalGndA,aNormalGndB))*180/M_PI); // rejected for normal
                    #endif
                    continue;
                }
                #ifdef SCANSCANDEBUG
                std::cout<<"Patch "<<aPatch.mLPatchesP[0]<<" accepted. Normals: "<<aNormalGndA<<" "<<aNormalGndB
                         <<" "<<acos(Scal(aNormalGndA,aNormalGndB))*180/M_PI<<"deg\n";
                #endif
                #ifdef SCANSCANSHOWPATCHES
                if (mBA.Iter()==mBA.NbMaxIter()-1)
                    mMapPatchesRasters.at(aScanA->NameImage()+"_to_"+aScanB->NameImage()).DIm().SetV(
                        aPatch.mLPatchesP[0]/SCANSCANSHOWPATCHES, aResidual); // accepted, give residual
                #endif
                aAvgRes.Add(1.0,fabs(aResidual));  // compute std deviation
                aVData.push_back(aData); // memorize the data for this image
            }
        } else {
            //std::cout<<" not visible\n";
            #ifdef SCANSCANSHOWPATCHES
            if (mBA.Iter()==mBA.NbMaxIter()-1)
                mMapPatchesRasters.at(aScanA->NameImage()+"_to_"+aScanB->NameImage()).DIm().SetV(
                    aPatch.mLPatchesP[0]/SCANSCANSHOWPATCHES, -20); // not visible
            #endif
        }
    }

    // if less than 1 scan to: nothing valuable to do
    if (aVData.size()<1) return NAN;

    mNbUsedPoints++;
    mNbUsedObs+=aVData.size();

    // accumlulate for computing average of deviation
    // mLastResidual.Add(1.0,  (aStdDev.StdDev(1e-5) *aVData.size()) / (aVData.size()-1.0));
    // mLastResidual.Add(1.0,  (aStdDev.StdDev(1e-5) ) );
    mLastResidual.Add(aVData.size(),  Square(aAvgRes.Average() ) );

    AddPatchDist(aPGround,aVData);

    for (auto & aData: aVData)
    {
        auto aCplId = std::pair(aData.mScanAName,aData.mScanBName);
        if (mMapNbUsedPatches.count(aCplId)==0)
            mMapNbUsedPatches[aCplId] = 1;
        else
            mMapNbUsedPatches[aCplId]++;
    }

    return aMinResidual;
}


void cBA_LidarLidarRaster::AddPatchDist
    (const cPt3dr & aPGround,
     const std::vector<cData1ImLidPhgr> &aVData
     )
{
    // read the solver now, because was not initialized at creation
    cResolSysNonLinear<tREAL8> *  aSys = mBA.Sys();
    // parse the data of the patch
    for (const auto & aData : aVData)
    {
        std::vector<int>       aVIndUk;
        std::vector<tREAL8>    aVObs;
        SetVUkVObs (aPGround,&aVIndUk,aVObs,aData,0);
        aSys->CalcAndAddObs(mEq,aVIndUk,aVObs,
                            *mWeightersMap.at(aData.mScanAName+"-"+aData.mScanBName));
    }
}



};
