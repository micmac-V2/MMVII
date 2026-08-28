#include "BundleAdjustment.h"
#include "../src/Graphs/ArboTriplets.h"

#include <latch>
#include <barrier>
#include <thread>
#include <mutex>
#include <condition_variable>

struct Barrier {
    int n, count = 0, generation = 0;
    std::mutex m;
    std::condition_variable cv;

    explicit Barrier(int n) : n(n) {}

    void wait() {
        std::unique_lock lock(m);
        int gen = generation;
        if (++count == n) {
            ++generation;
            count = 0;
            cv.notify_all();
        } else {
            cv.wait(lock, [&]{ return generation != gen; });
        }
    }
};



namespace MMVII
{

/*
class cThreadCreateSim
{
    public :
     
    private :
};
*/

Barrier aBar(5);

void MThreadForcConcurenceLock()
{
  if (!UserIsMPD()) return;
 // static int aCpt=0;

  aBar.wait();

  std::cout << "LTTTT " << TreeThreadsBase::Id() << " MPD?=" << UserIsMPD() << "\n";
}




/* ********************************************************* */
/*                                                           */
/*                     cBA_ArboTriplets                      */
/*                                                           */
/* ********************************************************* */

cBA_ArboTriplets::cBA_ArboTriplets(cMakeArboTriplet* aPMAT, std::vector<cSolLocNode>& aLocSols, int aTDepth, int aNbIterEnd,
                     tREAL8 aSigLooseningMult, tREAL8 aThrLooseningMult,cComputeMergeMulTieP * aTPts) :
    mPMAT      (aPMAT),
    mNbIter    (aNbIterEnd),
    mSigARange  ({2*aSigLooseningMult*aPMAT->Cfg().mSigmaAtt, Sqrt(aSigLooseningMult)*aPMAT->Cfg().mSigmaAtt}), //{max,min} <=> {initial,final}
    mThrRange   ({2*aThrLooseningMult*aPMAT->Cfg().mThrs,     aThrLooseningMult*aPMAT->Cfg().mThrs}),   //{max,min} <=> {initial,final}
    mSys      (nullptr),
    mTPts     (nullptr),
    mOwnTPts  (false),
    mTreeDepth(aTDepth)
{
    // get image names in current node
    std::vector<std::string> aVNames;
    for (auto & aSol : aLocSols)
        aVNames.push_back(aPMAT->MapI2Str(aSol.mNumPose));
    Sort2VectFirstOne(aVNames, aLocSols);

    // recover tie-points corresponding to the set of images
    if (aTPts)
    {
        // built once per node in MergeChildrenSol and shared with the similitude estimation
        mTPts = aTPts;
        MMVII_INTERNAL_ASSERT_medium(mTPts->VNames()==aVNames,
                                     "cBA_ArboTriplets : tie-points do not match the images of the node");
    }
    else
    {
        mTPts = new cComputeMergeMulTieP(*mPMAT->TPtsStruct(), aVNames);
        mOwnTPts = true;
    }

    // push initial values of intrinsics for your image set
    for (auto & aSol : aLocSols)
    {
        std::string aImName = aPMAT->MapI2Str(aSol.mNumPose);
        cPerspCamIntrCalib* aCal = aPMAT->PhProj().InternalCalibFromStdName(aImName, false);
        cIsometry3D<tREAL8> aPose(aSol.mPose.Tr(), aSol.mPose.Rot());
        mVCams.push_back(new cSensorCamPC(aImName, aPose, aCal));
        mVSens.push_back(mVCams.back());
        mVEqCol.push_back(mVCams.back()->CreateEqColinearityOnBundle(true, 100, false));
        mSetIntervUK.AddOneObj(mVCams.back());

    }
    // set-up the least squares system of equation
    mSys = new cResolSysNonLinear<tREAL8>(eModeSSR::eSSR_LsqNormSparse, mSetIntervUK.GetVUnKnowns());
    // vector of bundles decomposed to orthogonal u,v vectors
    mVecConfUV.resize(mTPts->Pts().size());
}

void cBA_ArboTriplets::OneIteration(int aIter)
{
    // viscosity
    for (auto& aCam : mVCams)
    {
        if ( mPMAT->Cfg().mViscPose.at(0)>0)
        {
            mSys->AddEqFixCurVar(*aCam,aCam->Center(),Square(1.0/mPMAT->Cfg().mViscPose.at(0)));
        }
        if (mPMAT->Cfg().mViscPose.at(1)>0)
        {
            mSys->AddEqFixCurVar(*aCam,aCam->Omega(),Square(1.0/mPMAT->Cfg().mViscPose.at(1)));
        }
    }

    // 3D intersection
    for (auto& aPair : mTPts->Pts())
        MakePGroundFromBundles(aPair, mVSens);

    // adapt the BA weight parameters
    if (aIter==0)
        AdaptWeightingToData();

    // diagnostic: compare triangulated P3D with GT at first iteration
    if (aIter==0 && mGTPts3D)
    {
        double aTotDist=0; int aNComp=0;
        double aMaxDist=0;
        for (auto& aAllConfigs : mTPts->Pts())
        {
            const auto& aVals = aAllConfigs.second;
            size_t aNbPts = aVals.mVIdPts.empty() ? NbPtsMul(aAllConfigs) : aVals.mVIdPts.size();
            for (size_t aKPts=0; aKPts<aNbPts; aKPts++)
            {
                if (aVals.mVIdPts.empty()) continue;
                int anId = aVals.mVIdPts.at(aKPts);
                auto it = mGTPts3D->find(anId);
                if (it == mGTPts3D->end()) continue;
                tREAL8 aDist = Norm2(aVals.mVPGround.at(aKPts) - it->second);
                aTotDist += aDist; aNComp++;
                UpdateMax(aMaxDist, (double)aDist);
            }
        }
        if (aNComp>0)
            StdOut() << "[DiagP3D] GT-vs-triangulated: avg=" << aTotDist/aNComp
                     << " max=" << aMaxDist << " over " << aNComp << " pts\n";
        else
            StdOut() << "[DiagP3D] no matching GT pts found (mVIdPts empty or no overlap)\n";
    }

    auto CurrentVal = [&](int aIterCur,int aIterMax,tREAL8 aV0,tREAL8 aV1)
    {
        if (aIterMax<=1) return aV1;
        return aV0 * std::pow(aV1/aV0, double(aIterCur)/(aIterMax-1));
    };

    tREAL8 aSigA = CurrentVal(aIter,mNbIter,mSigARange.at(0),mSigARange.at(1));
    tREAL8 aThr = CurrentVal(aIter,mNbIter,mThrRange.at(0),mThrRange.at(1));
    cStdWeighterResidual aTPtsW(1.0, aSigA, aThr, 2.0);

    //StdOutLock::lock();
    //StdOut() << "SIIGGGG THRRR " << aSigA << " " << aThr << std::endl;
    //StdOutLock::unlock();

    // add observation equations for all tie-points
    tREAL8 aMaxRes=0;
    int aNumAllTiePts=0;
    int aNumTPts=0;
    int aNumElimDegVis=0;   // eliminated by DegreeVisibility <= 0
    int aNumElimWeight=0;   // eliminated by weight == 0 (DegreeVisibility was > 0)
    cWeightAv<tREAL8> aWeigthedRes;

    int aConfigNum=0; //track id of current config

    // for every configuration of tie-pts
    for (auto aAllConfigs : mTPts->Pts())
    {
        const auto & aConfig = aAllConfigs.first;
        auto & aVals = aAllConfigs.second;

        size_t aNbIm = aConfig.size();
        // mVIdPts is only filled when created from MulTieP (with index); fall back to geometry-based count
        size_t aNbPts = aVals.mVIdPts.empty() ? NbPtsMul(aAllConfigs) : aVals.mVIdPts.size();


        // for every tie-point in current config
        for (size_t aKPts=0; aKPts<aNbPts; aKPts++)
        {

            const cPt3dr & aP3D = aVals.mVPGround.at(aKPts);
            cSetIORSNL_SameTmp<tREAL8>  aStrSubst(aP3D.ToStdVector());


            size_t aNbEqAdded = 0;

            // for every image where current tie-points is observed
            for (size_t aKIm=0; aKIm<aNbIm; aKIm++)
            {
                size_t aKImSorted = aConfig.at(aKIm);

                const cPt3dr aPBun (aVals.mVPIm.at(aKPts*aNbIm+aKIm).x(),
                                   aVals.mVPIm.at(aKPts*aNbIm+aKIm).y(),
                                   aVals.mVPZ.at(aKPts*aNbIm+aKIm) ) ;
                cSensorCamPC* aCam = this->mVCams.at(aKImSorted);

                // compute u,v orthogonal to bundle
                //    (do it once at first iteration)
                if (aIter==0)
                {
                    cPt3dr u = VUnit(VOrthog(aPBun));
                    cPt3dr v = VUnit(aPBun ^ u);
                    this->mVecConfUV.at(aConfigNum).push_back(std::make_pair(u,v));
                }

                // handle visibility
                //
              //   MThreadForcConcurenceLock();
                tREAL8 aDegVis = aCam->DegreeVisibility(aP3D);
                if (aDegVis > 0.0)
                {
                    // angle-based residual
                    const auto& [u, v] = mVecConfUV.at(aConfigNum).at(aKPts*aNbIm+aKIm);
                    cPt3dr aPBunPredUnit = VUnit(aCam->Pt_W2L(aP3D));

                    tREAL8 aF = aCam->InternalCalib()->F();
                    cPt2dr aResidual { aF * aPBunPredUnit.ToVect().DotProduct(u.ToVect()),
                                       aF * aPBunPredUnit.ToVect().DotProduct(v.ToVect()) };
                    tREAL8 aResNorm = Norm2(aResidual);

                    tREAL8 aWeight = aTPtsW.SingleWOfResidual(aResidual);

                    /*if (aWeight<=0.0)
                        StdOut() << "RRRR " << aResidual << " W=" << aWeight << ", "
                             << ", 3D=" << aP3D <<  ", F="
                             << aCam->InternalCalib()->F() << " "
                             << aPBun << " "
                             << aPBunPredUnit << "\n";*/


                    cCalculator<double> * aEqCol =  this->mVEqCol.at(aKImSorted);//aKIm


                    // add observations:
                    //    u,v vectors and focal
                    //    (rot init added implicitly in PushOwnObsColinearity)
                    std::vector<double> aVObs;
                    aVObs.push_back(mVecConfUV.at(aConfigNum).at(aKPts*aNbIm+aKIm).first.x());  //ux
                    aVObs.push_back(mVecConfUV.at(aConfigNum).at(aKPts*aNbIm+aKIm).first.y());  //uy
                    aVObs.push_back(mVecConfUV.at(aConfigNum).at(aKPts*aNbIm+aKIm).first.z());  //uz
                    aVObs.push_back(mVecConfUV.at(aConfigNum).at(aKPts*aNbIm+aKIm).second.x()); //vx
                    aVObs.push_back(mVecConfUV.at(aConfigNum).at(aKPts*aNbIm+aKIm).second.y()); //vy
                    aVObs.push_back(mVecConfUV.at(aConfigNum).at(aKPts*aNbIm+aKIm).second.z()); //vz
                    aVObs.push_back(aCam->InternalCalib()->F()); // focal

                    aCam->PushOwnObsColinearity(aVObs,aP3D);

                    std::vector<int> aVIndGlob = {-1,-2,-3};  // index of unknown, temporary
                    for (auto & anObj : aCam->GetAllUKPose())  // now add unknowns for sensor's extrinsics (no calib)
                    {
                        anObj->PushIndexes(aVIndGlob);
                    }

                    if (aWeight>0)
                    {
                        aWeigthedRes.Add(aWeight,aResNorm);//
                        mSys->R_AddEq2Subst(aStrSubst,aEqCol,aVIndGlob,aVObs,aWeight);//
                        aNbEqAdded++;
                        aNumTPts++;

                        if (aMaxRes<aResNorm)
                            aMaxRes=aResNorm;
                    }
                    else
                        aNumElimWeight++;
                }
                else
                    aNumElimDegVis++;
                aNumAllTiePts++;
            }

            if (aNbEqAdded>=2)
            {
                mSys->R_AddObsWithTmpUK(aStrSubst,mPMAT->Cfg().mLVM);
            }

        }
        aConfigNum++;
    }

    // print BA message only if at the last iterations
    if (aIter==(mNbIter-1))
    {
        StdOutLock::lock();
        StdOut() << "----------------------   Tree depth=" << mTreeDepth << ", images "
                 << NbCams() << "/" << mPMAT->GOP().AllVertices().size() << std::endl;

        StdOut() << "  # End Iter" << aIter+1
                 << " : Weighted Residual=" << aWeigthedRes.Average()
                 << ", #tie-points=" << aNumTPts << " #eliminated="
                 << " " << (100.0*(aNumAllTiePts-aNumTPts)/std::max(1,aNumAllTiePts)) << "%"
                 << " [DegVis<=0: " << aNumElimDegVis << ", Weight==0: " << aNumElimWeight << "]"
                 << std::endl;
        StdOutLock::unlock();
    }

    const auto& aVectSol = mSys->SolveUpdateReset({mPMAT->Cfg().mLVM}, {}, {});
    mSetIntervUK.SetVUnKnowns(aVectSol);
}

void cBA_ArboTriplets::UpdateLocSols(std::vector<cSolLocNode>& aLocSols)
{
    for (size_t aK = 0; aK < mVCams.size(); aK++)
    {
        aLocSols.at(aK).mPose.Tr()  = mVCams.at(aK)->Center();
        aLocSols.at(aK).mPose.Rot() = mVCams.at(aK)->Pose().Rot();
    }
}

tREAL8 cBA_ArboTriplets::RobustResidualScale(size_t aNbSample)
{
    size_t aNbObs = 0;
    for (const auto &aAllConfigs : mTPts->Pts())
        aNbObs+=aAllConfigs.second.mVPIm.size();

    const size_t aStep = std::max<size_t>(1,(aNbObs+aNbSample-1)/aNbSample);

    std::vector<double> aVRes;
    aVRes.reserve(std::min(aNbObs,aNbSample));
    size_t aKObs=0;

    for(const auto &aAllConfigs : mTPts->Pts())
    {
        const auto &aConf = aAllConfigs.first;
        const auto &aVals = aAllConfigs.second;
        const size_t aNbIm = aConf.size();
        const size_t aNbPts = aVals.mVPGround.size();

        for (size_t aKPts=0; aKPts<aNbPts; aKPts++)
        {
            const cPt3dr& aP3D = aVals.mVPGround.at(aKPts);
            for (size_t aKIm=0; aKIm<aNbIm; aKIm++,aKObs++)
            {
                if (aKObs % aStep) continue; /// take only N samples

                cSensorCamPC *aCam = mVCams.at(aConf.at(aKIm));
                if (aCam->DegreeVisibility(aP3D)<=0) continue; /// point must be visible

                const cPt3dr aPBun( aVals.mVPIm.at(aKPts*aNbIm+aKIm).x(),
                                   aVals.mVPIm.at(aKPts*aNbIm+aKIm).y(),
                                   aVals.mVPZ .at(aKPts*aNbIm+aKIm) );

                // |Unit(Bundle_obs) ^ Unit(Bundle_pred)| = sin(angle): exactly the norm
                // of the (u,v) residual of OneIteration, without needing u,v
                tREAL8 aSinA = Norm2( VUnit(aPBun) ^ VUnit(aCam->Pt_W2L(aP3D)) );
                aVRes.push_back(aCam->InternalCalib()->F() * aSinA);

            }
        }
    }
    if (aVRes.size() < 50) return -1; // not enough data => keep nominal weighting
    return NC_KthVal(aVRes,0.75);      // quantile, in pixels
}

void cBA_ArboTriplets::AdaptWeightingToData()
{
    const tREAL8 aSigAtt = mPMAT->Cfg().mSigmaAtt;
    const tREAL8 aThr    = mPMAT->Cfg().mThrs;

    /// compute quantile residual over this merged node
    mResScale = RobustResidualScale();
    if (mResScale<=0) return;

    static constexpr tREAL8 TheMaxLoosening = 40.0; /// not too sure about this value
    const tREAL8 aMult = std::clamp(mResScale/aSigAtt,1.0,TheMaxLoosening);

    StdOutLock::lock();
    StdOut() << aMult << " " << mResScale << " " << aSigAtt << std::endl;
    StdOutLock::unlock();

    /// update the SigmaAtt and Threshold ranges
    mSigARange = {2.0*aMult*aSigAtt, aSigAtt};
    mThrRange = {2*aMult*aThr, aThr};
}

cBA_ArboTriplets::~cBA_ArboTriplets()
{
    mSetIntervUK.SIUK_Reset();
    delete mSys;

    for (auto p : mVEqCol) delete p;
    for (auto p : mVCams)  delete p;

    if (mOwnTPts) delete mTPts;
}

};
