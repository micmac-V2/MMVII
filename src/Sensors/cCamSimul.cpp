#include "MMVII_PoseRel.h"
#include "MMVII_Tpl_Images.h"
#include "../Graphs/ArboTriplets.h"
#include "../BundleAdjustment/BundleAdjustment.h"

#include <unordered_set>

/* We have the image formula w/o distorsion:


*/

namespace MMVII
{

/* ************************************** */
/*                                        */
/*            cCamSimul                   */
/*                                        */
/* ************************************** */


cCamSimul::cCamSimul() :
   mCenterGround (10.0,5.0,20.0),
   mProfMin      (10.0),
   mProfMax      (20.0),
   mBsHMin       (0.1),
   mBsHMax       (0.5),
   mRandInterK   (0.1)
{
}

cCamSimul::~cCamSimul()
{
    DeleteAllAndClear(mListCam);
    DeleteAllAndClear(mListCalib);
}

const std::vector<cSensorCamPC *> & cCamSimul::ListCam() const { return mListCam; }

bool cCamSimul::ValidateCenter(const cPt3dr & aP2) const
{
    if (mListCam.empty()) return true;

    tREAL8 aTetaMin = 1e10;
    cPt3dr aV20 = aP2 - mCenterGround;
    for (const auto & aPtr : mListCam)
    {
         cPt3dr aV10 = aPtr->Center() - mCenterGround;
         UpdateMin(aTetaMin,AbsAngleTrnk(aV10,aV20));
    }
    return  (aTetaMin>mBsHMin) && (aTetaMin<mBsHMax);
}

cPt3dr  cCamSimul::GenCenterWOCstr(bool SubVert) const
{
    // Case "sub-vertical" we generate a point above mCenterGround
    //   * the delta in x and y is in an interval {
    if (SubVert)
    {
        auto v1 = RandUnif_C();
        auto v2 = RandUnif_C();
        auto v3 = RandInInterval(mProfMin,mProfMax);
        return    mCenterGround  + cPt3dr(v1,v2,1.0) * v3;
    }

    //
    auto v1 = cPt3dr::PRandUnit();
    auto v2 = RandInInterval(mProfMin,mProfMax);
    return mCenterGround + v1 * v2;
}


cPt3dr  cCamSimul::GenValideCenter(bool SubVert) const
{
   cPt3dr aRes = GenCenterWOCstr(SubVert);
   while (! ValidateCenter(aRes))
          aRes = GenCenterWOCstr(SubVert);

  // MMVII_INTERNAL_ASSERT_strong(!SubVert,"GenValideCenter");
  // StdOut() << "GenValideCenterGenValideCenter " << SubVert << "\n";
   return aRes;
}


void cCamSimul::AddCam(cPerspCamIntrCalib * aPC,bool SubVert,std::string Name)
{
      cPt3dr aNewC = GenValideCenter(SubVert);

      //  Axe K will point to the center of the scene
      cPt3dr aK = VUnit(mCenterGround - aNewC);
      // generate random I orthog to K
      cPt3dr aI = cPt3dr::PRandUnitNonAligned(aK,1e-2);
      // complete
      cPt3dr aJ = VUnit(aK ^aI);
      aI = aJ ^aK; // Just in case

      // we have now a rotation
      cRotation3D<tREAL8> aRot(M3x3FromCol(aI,aJ,aK),false);

      // if we add a small noise to not have a perfect intersec
      aNewC += cPt3dr::PRandC() * mRandInterK;
      // now we have a pose
      cIsometry3D<tREAL8> aPose(aNewC,aRot);

      // now we have a Cam
      mListCam.push_back(new cSensorCamPC(Name,aPose,aPC));
}

void cCamSimul::AddCam(eProjPC aProj,bool SubVert,std::string Name)
{
    // 1 => means Deg of direct dist is 2 (dir inverse is 5,1,1)
    cPerspCamIntrCalib * aCalib = cPerspCamIntrCalib::RandomCalib(aProj,1);

    mListCalib.push_back(aCalib);
    AddCam(aCalib,SubVert,Name);
}

cCamSimul * cCamSimul::Alloc2VIewTerrestrial(eProjPC aProj1,eProjPC aProj2,bool SubVert)
{
   cCamSimul * aRes = new cCamSimul();

   aRes->AddCam(aProj1,SubVert);
   aRes->AddCam(aProj2,SubVert);

   return aRes;
}

cCamSimul * cCamSimul::AllocNVIewTerrestrial(int aNb,eProjPC aProj,bool SubVert)
{
    cCamSimul * aRes = new cCamSimul();

    for (int aK=0; aK<aNb; aK++)
    {
        aRes->AddCam(aProj,SubVert,"SimCam"+ToStr(aK));
    }

    return aRes;
}

void cCamSimul::TestCam(cSensorCamPC * aCam) const
{
        StdOut() << "CC " << aCam->Center()  << " CG=" << mCenterGround << std::endl;

cPt3dr aV = aCam->Center() - mCenterGround;

StdOut()  << " I " << Cos(aV,aCam->AxeI())
          << " J " << Cos(aV,aCam->AxeI())
          << " K " << Cos(aV,aCam->AxeK())
          << " \n";

        StdOut() << "Vis " <<  aCam->IsVisible(mCenterGround) << std::endl;
}

void BenchMEP_Coplan();

void cCamSimul::BenchPoseRel2Cam
     (
        cTimerSegm * aTS,
        bool         PerfInter,
        bool         isSubVert,
        bool         isPlanar
     )
{
    thread_local static int aCpt=0;
    /// cLinearOverCstrSys<tREAL8> *  aSysL1 = AllocL1_Barrodale<tREAL8>(9);
    // cLinearOverCstrSys<tREAL8> *  aSysL1 = new cLeasSqtAA<tREAL8>(9);
    cLeasSqtAA<tREAL8> aSysL2(9);

    thread_local static int aCptPbL1 = 0; FakeUseIt(aCptPbL1);


    if (1)
    {

       for (int aK1=0 ; aK1<(int)eProjPC::eNbVals ; aK1++)
       {
           for (int aK2=0 ; aK2<(int)eProjPC::eNbVals ; aK2++)
           {
               cAutoTimerSegm aTSSim(aTS,"CreateSimul");
               aCpt++;
               cCamSimul * aCamSim = cCamSimul::Alloc2VIewTerrestrial(eProjPC(aK1),eProjPC(aK2),isSubVert);

               // we want to test robustness in perfect degenerate & close to degenertae
               if (PerfInter)
                  aCamSim->mRandInterK = 0.0;

               // Generate 2 cams
               cSensorCamPC * aCam1 = aCamSim->mListCam.at(0);
               cSensorCamPC * aCam2 = aCamSim->mListCam.at(1);

               // generate  perfect homologous point
               cSetHomogCpleIm aSetH;
               size_t aNbPts = 40;

               for (size_t aKP=0 ; aKP<aNbPts ; aKP++)
               {
                  // StdOut() << " Planaaarr " << isPlanar << " K=" << aKP << "\n";
                  cHomogCpleIm aCple =  isPlanar                                                     ?
                                        aCam1->RandomVisibleCple(aCamSim->mCenterGround.z(),*aCam2)  :
                                        aCam1->RandomVisibleCple(*aCam2)                             ;
                  aSetH.Add(aCple);
               }

      ///     StdOut() << "Ouut IssssPllannnn " << isPlanar << "\n";

               // Make 3D direction of points
               cSetHomogCpleDir aSetD (aSetH,*(aCam1->InternalCalib()),*(aCam2->InternalCalib()));

               cAutoTimerSegm aTSGetMax(aTS,"GetMaxK");
         //   StdOut() << "Ouut IssssPllannnn " << isPlanar << "\n";

               if (isPlanar )
               {
                       // To see, not sure validate any more
                       /*
                    cPSC_PB aParam("Cam");
                    cPS_CompPose aPsC(aSetD,&aParam);
                    */
               }
               else if (false)
               {
                   int aKMax =  MatEss_GetKMax(aSetD,1e-6);

                  // These point where axe k almost intersect, the z1z2 term of mat ess is probably small
                  // and must not be KMax
                   MMVII_INTERNAL_ASSERT_bench(aKMax!=8,"cComputeMatEssential::GetKMax");

                // Now test that residual is ~ 0 on these perfect points
                   cAutoTimerSegm aTSL2(aTS,"L2");
                   cMatEssential aMatEL2(aSetD,aSysL2,aKMax);

                   {
                       cIsometry3D<tREAL8>  aPRel =  aCam1->RelativePose(*aCam2);
                       // When we give aPRel
                       aMatEL2.ComputePose(aSetD,&aPRel);
                   }
                   MMVII_INTERNAL_ASSERT_bench(aMatEL2.AvgCost(aSetD,1.0)<1e-5,"Avg cost ");

                  cAutoTimerSegm aTSL1(aTS,"L1");
                  cLinearOverCstrSys<tREAL8> *  aSysL1 = AllocL1_Barrodale<tREAL8>(9);
                  cMatEssential aMatEL1(aSetD,*aSysL1,aKMax);
                  MMVII_INTERNAL_ASSERT_bench(aMatEL1.AvgCost(aSetD,1.0)<1e-5,"Avg cost ");

                  for (int aK=0 ; aK<4 ; aK++)
                      aSetD.GenerateRandomOutLayer(0.1);

                  cMatEssential aMatNoise(aSetD,*aSysL1,aKMax);

                  delete aSysL1;
            
                 if (0)
                 {
                     StdOut() << "Cpt=" << aCpt
                         << " Cost95= "  << aMatNoise.KthCost(aSetD,0.95)
                         << " Cost80= "  << aMatNoise.KthCost(aSetD,0.70)
                         << " KMax= "  << aKMax
                         << "\n" ;
                     MMVII_INTERNAL_ASSERT_bench(aMatNoise.KthCost(aSetD,0.70) <1e-5,"Kth cost ");
                     MMVII_INTERNAL_ASSERT_bench(aMatNoise.KthCost(aSetD,0.95) >1e-2,"Kth cost ");
                 }

                 // We test if the residual at 70% is almost 0 (with 4/40 outlayers)
                 if (aMatNoise.KthCost(aSetD,0.70)>1e-5)
                      aCptPbL1++;
               }
               delete aCamSim;
            }
        }
    }

    BenchMEP_Coplan();
}


struct TripletHash {
    std::size_t operator()(const std::array<int,3>& t) const noexcept {
        std::size_t h1 = std::hash<int>{}(t[0]);
        std::size_t h2 = std::hash<int>{}(t[1]);
        std::size_t h3 = std::hash<int>{}(t[2]);

        // hash combine
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

void Bench_HBA(cParamExeBench & aParam)
{
    if (! aParam.NewBench("HierarchBA")) return;

    cTimerSegm * aTS = nullptr;
    if (aParam.Show())
    {
        aTS = new cTimerSegm(&cMMVII_Appli::CurrentAppli());
    }

    cCamSimul::BenchHierchBA_InitOnly(aTS,false);
    cCamSimul::BenchHierchBA_BAOnly(aTS,false);
    //cCamSimul::BenchHierchBA(aTS,true,false);

    delete aTS;
    aParam.EndBench();
}

/*  Shared helpers for the dissected benchmarks   */

struct cBenchScene
{
    std::unique_ptr<cCamSimul>      mCamSim;
    std::vector<cDataSolOriTriplet> m3Set;
    std::vector<std::string>        mSetIm;
    std::map<int,cPt3dr>            mGTPts3D;  ///< GT 3D coords keyed by tie-point ID (aPtImIdx)
};

/// Generate cameras, triplets, tie-points and GT relative poses for one projection type
static cBenchScene BuildBenchScene(int aNbCam, int aNbTri, int aNbHPts,
                                   eProjPC aProj, bool isSubVert,
                                   double aPtNoiseAmpl,
                                   double aPtOutlierRate, double aPtOutlierAmpl,
                                   double aTriOutNb, const std::vector<double>& aTriOutAmpl,
                                   cPhotogrammetricProjectMemory& aMemPhProj)
{
    cBenchScene aS;
    aS.mCamSim.reset(cCamSimul::AllocNVIewTerrestrial(aNbCam, aProj, isSubVert));
    // mRandInterK is private - caller (a static member) sets it after construction

    const auto& aListCam = aS.mCamSim->ListCam();
    for (int aKIm = 0; aKIm < aNbCam; aKIm++)
        aMemPhProj.AddCalib(aListCam.at(aKIm)->NameImage(),
                            aListCam.at(aKIm)->InternalCalib());

    std::vector<std::pair<int,int>> aEdges;
    std::vector<bool> aCamVisited(aNbCam, false);
    std::unordered_set<std::array<int,3>, TripletHash> aTriplets;

    while (aTriplets.empty())
    {
        auto aFirstTri = RandSet(3, aNbCam);
        std::array<int,3> t = {aFirstTri[0], aFirstTri[1], aFirstTri[2]};
        std::sort(t.begin(), t.end());
        if (t[0]==t[1] || t[0]==t[2] || t[1]==t[2]) continue;
        aTriplets.insert(t);
        for (int n : t) { aCamVisited[n]=1; aEdges.push_back({t[0],t[1]}); aEdges.push_back({t[1],t[2]}); aEdges.push_back({t[2],t[0]}); }
    }
    while ((int)aTriplets.size() < aNbTri)
    {
        int anEdgeId = RandUnif_N(aEdges.size()-1);
        std::array<int,3> t = {aEdges[anEdgeId].first, aEdges[anEdgeId].second, (int)RandUnif_N(aNbCam-1)};
        std::sort(t.begin(), t.end());
        if (t[0]==t[1] || t[0]==t[2] || t[1]==t[2]) continue;
        auto [it, inserted] = aTriplets.insert(t);
        if (inserted)
            for (int n : t) { aCamVisited[n]=1; aEdges.push_back({t[0],t[1]}); aEdges.push_back({t[1],t[2]}); aEdges.push_back({t[2],t[0]}); }
    }

    for (size_t aKCV=0; aKCV<aCamVisited.size(); aKCV++)
        if (aCamVisited[aKCV]) aS.mSetIm.push_back(aListCam[aKCV]->NameImage());
    std::sort(aS.mSetIm.begin(), aS.mSetIm.end());

    std::map<std::string, cVecTiePMul> aMulTiePMap;
    for (auto* aCam : aListCam) aMulTiePMap[aCam->NameImage()] = cVecTiePMul();

    int aPtImIdx = 0;
    for (auto& aT : aTriplets)
    {
        std::vector<cSensorCamPC*> aCams = {aListCam.at(aT[0]),
                                             aListCam.at(aT[1]),
                                             aListCam.at(aT[2])};
        int aNbKeyPtsInTri = 0;
        for (size_t aK1Cam=0; aK1Cam<aCams.size(); aK1Cam++)
            for (size_t aK2Cam=aK1Cam+1; aK2Cam<aCams.size(); aK2Cam++)
                for (int aKPt=0; aKPt<aNbHPts; aKPt++)
                {
                    // Generate P3D on cam1's ray so it lies exactly on that ray.
                    // PInterBundle of two random pixels gives the midpoint of skew rays,
                    // which is NOT on either ray for non-perspective cameras (e.g. eEquiRect),
                    // causing non-zero BA residuals even with GT poses.
                    cPt3dr aPt3D;
                    cHomogCpleIm aHPair;
                    bool aGotPair = false;
                    for (int aTry=0; aTry<1000 && !aGotPair; aTry++)
                    {
                        aPt3D = aCams[aK1Cam]->RandomVisiblePGround(5.0, 40.0);
                        cPt2dr aPIm2 = aCams[aK2Cam]->Ground2Image(aPt3D);
                        if (aCams[aK2Cam]->IsVisible(aPt3D) && aCams[aK2Cam]->IsVisibleOnImFrame(aPIm2))
                        {
                            aHPair = cHomogCpleIm(aCams[aK1Cam]->Ground2Image(aPt3D), aPIm2);
                            aGotPair = true;
                        }
                    }
                    if (!aGotPair) continue;
                    aS.mGTPts3D[aPtImIdx] = aPt3D;
                    for (size_t aK3Cam=0; aK3Cam<aCams.size(); aK3Cam++)
                        if (aK3Cam!=aK1Cam && aK3Cam!=aK2Cam && aCams[aK3Cam]->IsVisible(aPt3D))
                        {
                            cPt2dr aPt = aCams[aK3Cam]->Ground2Image(aPt3D);
                            cPt2dr aDelta = cPt2dr::PRandC()*aPtNoiseAmpl;
                            if (aCams[aK3Cam]->IsVisibleOnImFrame(aPt+aDelta)) aPt += aDelta;
                            aMulTiePMap[aCams[aK3Cam]->NameImage()].mVecTPM.push_back(cTiePMul(aPt, aPtImIdx));
                            aNbKeyPtsInTri++;
                        }
                    cPt2dr aDelta1 = cPt2dr::PRandC()*aPtNoiseAmpl;
                    cPt2dr aDelta2 = cPt2dr::PRandC()*aPtNoiseAmpl;
                    if (aCams[aK1Cam]->IsVisibleOnImFrame(aHPair.mP1+aDelta1)) aHPair.mP1 += aDelta1;
                    if (aCams[aK2Cam]->IsVisibleOnImFrame(aHPair.mP2+aDelta2)) aHPair.mP2 += aDelta2;
                    aMulTiePMap[aCams[aK1Cam]->NameImage()].mVecTPM.push_back(cTiePMul(aHPair.mP1, aPtImIdx));
                    aMulTiePMap[aCams[aK2Cam]->NameImage()].mVecTPM.push_back(cTiePMul(aHPair.mP2, aPtImIdx));
                    aNbKeyPtsInTri += 2;
                    aPtImIdx++;
                }

        int aNbOutliers = (int)(aNbKeyPtsInTri * aPtOutlierRate);
        for (int aKOut=0; aKOut<aNbOutliers; aKOut++)
        {
            int aRandCam = RandUnif_N(3);
            int aNbKP = (int)aMulTiePMap[aCams[aRandCam]->NameImage()].mVecTPM.size();
            int aRandIdx = RandUnif_N(aNbKP);
            cPt2dr aPt = aMulTiePMap[aCams[aRandCam]->NameImage()].mVecTPM.at(aRandIdx).mPt;
            cPt2dr aDelta = cPt2dr::PRandC()*aPtOutlierAmpl;
            if (aCams[aRandCam]->IsVisibleOnImFrame(aPt+aDelta))
                aMulTiePMap[aCams[aRandCam]->NameImage()].mVecTPM.at(aRandIdx).mPt = aPt+aDelta;
        }

        cSensorCamPC* aCam1 = aCams[0];
        cSensorCamPC* aCam2 = aCams[1];
        cSensorCamPC* aCam3 = aCams[2];
        tPoseR aPose2to1 = aCam1->RelativePose(*aCam2);
        tPoseR aPose3to1 = aCam1->RelativePose(*aCam3);
        double aDist = Norm2(aPose2to1.Tr());
        aPose2to1.Tr() = aPose2to1.Tr() / aDist;
        aPose3to1.Tr() = aPose3to1.Tr() / aDist;
        cDataSolOriTriplet aTri;
        aTri.mVNames = {aCam1->NameImage(), aCam2->NameImage(), aCam3->NameImage()};
        aTri.mP01 = aPose2to1;  aTri.mP02 = aPose3to1;
        aS.m3Set.push_back(aTri);
    }

    for (int aKtri=0; aKtri<(int)aTriOutNb; aKtri++)
    {
        int aRandTri = RandUnif_N((int)aS.m3Set.size());
        for (int aK=1; aK<3; aK++)
        {
            tPoseR& aCurP = (aK==1) ? aS.m3Set[aRandTri].mP01 : aS.m3Set[aRandTri].mP02;
            aCurP.Tr() += cPt3dr::PRandInSphere() * aTriOutAmpl[0];
            aCurP.Rot() = aCurP.Rot() * tRotR::RandomSmallElem(aTriOutAmpl[1]);
        }
    }

    for (auto& aHStr : aMulTiePMap)
        aMemPhProj.AddMulTieP(aHStr.first, aHStr.second);

    return aS;
}

/// Computes similarity-aligned pose errors vs GT
static double CompareWithGT(const std::vector<cSensorCamPC*>& aGTCams,
                            std::map<std::string, cSensorCamPC*>& aSolCams)
{
    std::vector<tPoseR> aVPoseGT, aVPoseCalc;
    for (auto* aCamGT : aGTCams)
        if (aSolCams[aCamGT->NameImage()])
        {
            aVPoseGT.push_back(aCamGT->Pose());
            aVPoseCalc.push_back(aSolCams[aCamGT->NameImage()]->Pose());
        }
    auto [aRes, aSim] = EstimateSimTransfertFromPoses(aVPoseGT, aVPoseCalc);

    double aErrTr=0, aErrR=0;
    int aN=0;
    for (auto* aCamGT : aGTCams)
        if (aSolCams[aCamGT->NameImage()])
        {
            tPoseR aPoseCalcInGT = TransfoPose(aSim, aSolCams[aCamGT->NameImage()]->Pose());
            aErrTr += Norm2(aCamGT->Pose().Tr() - aPoseCalcInGT.Tr());
            aErrR  += aCamGT->Pose().Rot().Dist(aPoseCalcInGT.Rot());
            aN++;
        }
    StdOut() << "ErrTrAvg=" << aErrTr/aN << ", ErrRAvg=" << aErrR/aN << "\n";
    return aErrTr / aN;
}

/*  Bench 0: spanning tree initialisation followed by BA */

void cCamSimul::BenchHierchBA(cTimerSegm * aTS,
                              bool         PerfInter,
                              bool         isSubVert)
{
    const int aNbCam = 25;
    const int aNbTri = round_ni((aNbCam-1)/2 * 5);
    const int aNbHPts = 20;

    cMakeArboTripletCfg aCfg;
    aCfg.mLVM      = 0.1;
    aCfg.mNbIterBA = 5;
    aCfg.mSigmaTPt = 1;
    aCfg.mFacElim  = 10;

    StdOut() << "Nb of cams=" << aNbCam << ", nb of triplets=" << aNbTri << std::endl;

    cMMVII_Appli& anAp = cMMVII_Appli::CurrentAppli();

    for (size_t aK1=0; aK1<(size_t)eProjPC::eNbVals; aK1++)
    {
        StdOut() << "=== Bench Init+BA,  Projection=" << E2Str(eProjPC(aK1)) << " ===\n";

        cPhotogrammetricProjectMemory aMemPhProj;
        cBenchScene aScene = BuildBenchScene(aNbCam, aNbTri, aNbHPts, eProjPC(aK1), isSubVert,
                                             0.0, 0.0, 20.0, 0.0, {0.5, 0.1}, aMemPhProj);
        if (PerfInter)
            aScene.mCamSim->mRandInterK = 0.0;

        StdOut() << "Start Hierarchical SfM" << std::endl;
        cMakeArboTriplet aMk3(aScene.m3Set, false, 1.0, aMemPhProj, anAp, aCfg);
        aMk3.InitTPtsStruct("", aScene.mSetIm);
        aMk3.MakeGraphPose();
        aMk3.InitialiseCalibs();
        aMk3.DoPoseRef();
        aMk3.MakeCnxTriplet();
        aMk3.MakeWeightingGraphTriplet();
        aMk3.ComputeArbor();
        aMk3.SaveGlobSol();

        auto aSolCams = aMemPhProj.SensorMap();
        CompareWithGT(aScene.mCamSim->ListCam(), aSolCams);
    }
}




/*  Bench 1: spanning tree initialisation only (no BA) */

void cCamSimul::BenchHierchBA_InitOnly(cTimerSegm* aTS, bool isSubVert)
{
    const int aNbCam = 25;
    const int aNbTri = round_ni((aNbCam-1)/2 * 5);
    const int aNbHPts = 20;

    cMakeArboTripletCfg aCfg;
    aCfg.mLVM      = 0.1;
    aCfg.mNbIterBA = 0;   // spanning tree only — no BA refinement
    aCfg.mSigmaTPt = 1;
    aCfg.mFacElim  = 10;

    cMMVII_Appli& anAp = cMMVII_Appli::CurrentAppli();

    for (size_t aK1=0; aK1<(size_t)eProjPC::eNbVals; aK1++)
    {
        StdOut() << "=== BenchInitOnly  Projection=" << E2Str(eProjPC(aK1)) << " ===\n";

        cPhotogrammetricProjectMemory aMemPhProj;
        cBenchScene aScene = BuildBenchScene(aNbCam, aNbTri, aNbHPts, eProjPC(aK1), isSubVert,
                                             0.0, 0.0, 20.0, 0.1, {0.5,0.1}, aMemPhProj);
        aScene.mCamSim->mRandInterK = 0.0;

        cMakeArboTriplet aMk3(aScene.m3Set, false, 1.0, aMemPhProj, anAp, aCfg);
        aMk3.InitTPtsStruct("", aScene.mSetIm);
        aMk3.MakeGraphPose();
        aMk3.InitialiseCalibs();
        aMk3.DoPoseRef();
        aMk3.MakeCnxTriplet();
        aMk3.MakeWeightingGraphTriplet();
        aMk3.ComputeArbor();
        aMk3.SaveGlobSol();

        auto aSolCams = aMemPhProj.SensorMap();
        CompareWithGT(aScene.mCamSim->ListCam(), aSolCams);
    }
}

/*  Bench 2: BA refinement only - GT poses as initialisation           */

void cCamSimul::BenchHierchBA_BAOnly(cTimerSegm* aTS, bool isSubVert)
{
    const int aNbCam = 25;
    const int aNbTri = round_ni((aNbCam-1)/2 * 5);
    const int aNbHPts = 20;
    const int aNbIterBA = 5;

    cMakeArboTripletCfg aCfg;
    aCfg.mLVM      = 0.1;
    aCfg.mNbIterBA = aNbIterBA;
    aCfg.mSigmaTPt = 1;
    aCfg.mFacElim  = 10;
    //aCfg.mViscPose = {0.1,0.1};

    cMMVII_Appli& anAp = cMMVII_Appli::CurrentAppli();

    for (size_t aK1=0; aK1<(size_t)eProjPC::eNbVals; aK1++)
    {
        StdOut() << "=== BenchBAOnly  Projection=" << E2Str(eProjPC(aK1)) << " ===\n";

        cPhotogrammetricProjectMemory aMemPhProj;
        cBenchScene aScene = BuildBenchScene(aNbCam, aNbTri, aNbHPts, eProjPC(aK1), isSubVert,
                                             0.0, 0.0, 20.0, 0, {0.5,0.1}, aMemPhProj);
        aScene.mCamSim->mRandInterK = 0.0;

        cMakeArboTriplet aMk3(aScene.m3Set, false, 1.0, aMemPhProj, anAp, aCfg);
        aMk3.InitTPtsStruct("", aScene.mSetIm);
        aMk3.MakeGraphPose();
        aMk3.InitialiseCalibs();

        // build GT initial poses for all cameras in the pose graph
        std::vector<cSolLocNode> aLocSols;
        for (int aKP=0; aKP<(int)aMk3.GOP().NbVertex(); aKP++)
        {
            const std::string& aName = aMk3.MapI2Str(aKP);
            for (auto* aGTCam : aScene.mCamSim->ListCam())
                if (aGTCam->NameImage() == aName)
                {
                    aLocSols.push_back(cSolLocNode(aGTCam->Pose(), aKP));
                    break;
                }
        }

        // run BA from GT initial poses
        cBA_ArboTriplets aBA(&aMk3, aLocSols);
        aBA.SetGTPts3D(&aScene.mGTPts3D);
        for (int aIter=0; aIter<aNbIterBA; aIter++)
            aBA.OneIteration(aIter);
        aBA.UpdateLocSols(aLocSols);

        // save refined poses so we can use CompareWithGT
        for (auto& aSol : aLocSols)
        {
            const std::string& aName = aMk3.MapI2Str(aSol.mNumPose);
            cPerspCamIntrCalib* aCal = aMk3.PhProj().InternalCalibFromImage(aName);
            cSensorCamPC aCamPC(aName, aSol.mPose, aCal);
            aMemPhProj.SaveCamPC(aCamPC);
        }

        auto aSolCams = aMemPhProj.SensorMap();
        double aErrTr = CompareWithGT(aScene.mCamSim->ListCam(), aSolCams);
        MMVII_INTERNAL_ASSERT_bench(aErrTr < 1e-3, "BenchHierchBA_BAOnly: BA from GT should give near-zero error on perfect data");
    }
}

}; // MMVII




