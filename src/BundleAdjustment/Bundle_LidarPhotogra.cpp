#include "BundleAdjustment.h"
#include "MMVII_Interpolators.h"
#include "MMVII_Tpl_Images.h"

namespace MMVII
{

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



cBA_LidarPhotogra::~cBA_LidarPhotogra()
{
}

cBA_LidarPhotograTri::~cBA_LidarPhotograTri()
{
    if (mTri) delete mTri;
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
             continue; // do not project to other TLS
          // 1st test: zbuffer visibility
          //std::cout<<"Im "<<aScanBData.mScanName<<" patch "<<aPatch.mId<<" vis "<<aPatch.mImVisible.at(aScanBData.mScanName)<<"\n";
          if (aHiddenOnImage.count(aCam->NameImage())>0)
          {
              continue;
          }
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
                        if (!aPIm.IsValid())
                        {
                           continue;
                        }
                        if (aGenDIm.InsideInterpolator(*mInterp,aPIm,1.0))  // is it sufficiently inside
                        {
                            auto aVGr = aGenDIm.GetValueAndGradInterpol(*mInterp,aPIm); // extract pair Value/Grad of image
                            if ((aVGr.first==0)||(aVGr.first==255)) // refuse saturated pixels TODO improve criteria!
                                continue;
                            aData.mVGr.push_back(aVGr); // push it at end of stack
                        } else {
                        }
                   } else {
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

};
