#include "BundleAdjustment.h"
#include "../src/Graphs/ArboTriplets.h"

namespace MMVII
{

/* ********************************************************* */
/*                                                           */
/*                     cBA_ArboTriplets                      */
/*                                                           */
/* ********************************************************* */

cBA_ArboTriplets::cBA_ArboTriplets(cMakeArboTriplet* aPMAT, std::vector<cSolLocNode>& aLocSols):
    mPMAT      (aPMAT),
    mNbIter    (aPMAT->NbIterBA()),
    mSigAttFinal(1.0),
    mThrFinal   (aPMAT->FacElim()),
    //mSigARange  ({std::max(mSigAttFinal,std::min(5.0,aPMAT->SigmaTPt())),mSigAttFinal}), // {max,min} <=> {initial,final}
    //mThrRange   ({std::max(mThrFinal,std::min(30.0,aPMAT->SigmaTPt()*aPMAT->FacElim())),mThrFinal}), // {max,min} <=> {initial,final}
    mSigARange  ({std::max(mSigAttFinal,aPMAT->SigmaTPt()),mSigAttFinal}), // {max,min} <=> {initial,final}
    mThrRange   ({mThrFinal,mThrFinal}),            // {max,min} <=> {initial,final}
    mSys      (nullptr),
    mTPts     (nullptr)
{
    // get image names in current node
    std::vector<std::string> aVNames;
    for (auto & aSol : aLocSols)
        aVNames.push_back(aPMAT->MapI2Str(aSol.mNumPose));
    Sort2VectFirstOne(aVNames, aLocSols);

    // recover tie-points corresponding to the set of images
    mTPts = new cComputeMergeMulTieP(*mPMAT->TPtsStruct(), aVNames);

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
        if ( mPMAT->ViscPose().at(0)>0)
        {
            mSys->AddEqFixCurVar(*aCam,aCam->Center(),Square(1.0/mPMAT->ViscPose().at(0)));
        }
        if (mPMAT->ViscPose().at(1)>0)
        {
            mSys->AddEqFixCurVar(*aCam,aCam->Omega(),Square(1.0/mPMAT->ViscPose().at(1)));
        }
    }

    // 3D intersection
    for (auto& aPair : mTPts->Pts())
        MakePGroundFromBundles(aPair, mVSens);

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

    auto CurrentVal = [&](int iterCur,int iterMax,tREAL8 delta,tREAL8 bias)
    {
        return delta*(1 - double(iterCur)/(iterMax-1)) + bias;
    };

    tREAL8 aSigA = CurrentVal(aIter,mNbIter,mSigARange.at(0) - mSigARange.at(1),mSigARange.at(1));
    tREAL8 aThr = CurrentVal(aIter,mNbIter,mThrRange.at(0) - mThrRange.at(1),mThrRange.at(1));
    cStdWeighterResidual aTPtsW(1.0, aSigA, aThr, 2.0);

    // add observation equations for all tie-points
    tREAL8 aMaxRes=0;
    int aNumAllTiePts=0;
    int aNumTPts=0;
    int aNumAll3DPts=0;
    int aNum3DPts=0;
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

        aNumAll3DPts+=aNbPts;

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
                mSys->R_AddObsWithTmpUK(aStrSubst,mPMAT->LVM());
                aNum3DPts++;
            }

        }
        aConfigNum++;
    }

    double aPercInliersTP = (aNumAllTiePts>0) ? (aNumTPts*100)/aNumAllTiePts : 0.0;
    double aPercIn3DP = (aNumAll3DPts>0) ? (aNum3DPts*100)/aNumAll3DPts : 0.0;
    if (aIter==0)
        StdOut() << "[iter0] residuals: weightedAvg=" << aWeigthedRes.Average()
                 << " max=" << aMaxRes
                 << " #eliminated=" << (aNumAllTiePts-aNumTPts) << "/" << aNumAllTiePts
                 << " (" << (100.0*(aNumAllTiePts-aNumTPts)/std::max(1,aNumAllTiePts)) << "%)"
                 << " [DegVis<=0: " << aNumElimDegVis << ", Weight==0: " << aNumElimWeight << "]\n";
    StdOut() << "#Iter=" << aIter
             << " Res=" << aWeigthedRes.Average()
             << ", #3D points=" << aNumAll3DPts << ", " << aPercIn3DP << "%"
             << ", #2D features=" << aNumTPts << ", " << aPercInliersTP << "%"
             // << ", MaxRes=" << aMaxRes
             << std::endl;

    const auto& aVectSol = mSys->SolveUpdateReset({mPMAT->LVM()}, {}, {});
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

cBA_ArboTriplets::~cBA_ArboTriplets()
{
    mSetIntervUK.SIUK_Reset();
    delete mSys;
    delete mTPts;
    for (auto p : mVEqCol) delete p;
    for (auto p : mVCams)  delete p;
}

};
