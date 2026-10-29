/** \file cAppli_SelecTripletsOriRel.cpp
    \brief Selection of triplets by coverage of image and B/H, port of MMV1 Martini NO_GenTripl
*/

#include "MMVII_2Include_Serial_Tpl.h"
#include "MMVII_DeclareAllCmd.h"
#include "MMVII_PoseRel.h"
#include "MMVII_Sensor.h"
#include "MMVII_MeasuresIm.h"
#include "SelecTripletsCoverage.h"
#include <algorithm>

namespace MMVII
{

/// coordinates used for the grid of coverage : calibrated ray coordinates (x/z,y/z) ; assume z>0 (not true for wide fisheye?)
static cPt2dr GridCoordOfDir(const cPt3dr & aDir)
{
    return cPt2dr(aDir.x()/aDir.z(),aDir.y()/aDir.z());
}

/* ********************************************************** */
/*                                                            */
/*                 cAppli_SelecTripletsOriRel                 */
/*                                                            */
/* ********************************************************** */

/** Selection of triplets for OriPoseEstimRelAllTriplets : for each oriented pair (A,B), A<B, greedy choice
 *  of third images C maximising the coverage of image A by triple tie points, weighted by B/H.
 *  Mode 1 : all oriented pairs (A,B) stored in folder of A ; mode 2 : all images (driver).
 */
class cAppli_SelecTripletsOriRel : public cMMVII_Appli
{
public :
    typedef std::pair<std::string,std::string> tPairNames;

    cAppli_SelecTripletsOriRel(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec,int aMode);
    int Exe() override;
    cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
    cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;

private :
    void ReadOrientedPairs();      ///< fill mVNamesIm, mMapPoses, mNeighbours
    void ReadCycleTriplets();      ///< fill mThirdsOfCycle, from the copy "Cyc." if it exists
    bool IsOriented(const std::string & aN1,const std::string & aN2) const;
    /// pose Im2->Im1 (first candidate), inverted if the pair is stored as (Im2,Im1)
    tPoseR PoseOfPair(const std::string & aN1,const std::string & aN2) const;
    /// tie points of a stored pair (aN1<aN2), filtered on residual, cached
    const cSetHomogCpleIm & FilteredHom(const std::string & aN1,const std::string & aN2);
    std::vector<std::string> Candidates(const std::string & aNA,const std::string & aNB) const;
    /// selected triplets of the oriented pair (aNA,aNB), aNA<aNB
    std::vector<cTripletName> SelecOfPair(const std::string & aNA,const std::string & aNB);

    int                       mModeCompute;
    cPhotogrammetricProject   mPhProj;
    std::string               mIm1;
    std::string               mIm2;           ///< mode 1, optional : process only pair (Im1,Im2), debug
    cParamSelTriCov           mParamSel;
    eTriCandPool              mCandPool;
    int                       mNbMax;
    int                       mMinNb3;
    tREAL8                    mPropResFilter;
    bool                      mShow;

    std::vector<std::string>                         mVNamesIm;       ///< all images of the OriRel folder
    std::map<tPairNames,std::vector<tPoseR>>         mMapPoses;       ///< oriented pairs (Im1<Im2), all candidate poses Im2->Im1
    std::map<std::string,std::set<std::string>>      mNeighbours;     ///< graph of oriented pairs
    std::map<tPairNames,std::set<std::string>>       mThirdsOfCycle;  ///< pool (ii) : pair -> third images of cycle triplets
    std::map<tPairNames,cSetHomogCpleIm>             mCacheHom;       ///< filtered tie points by stored pair
};

cAppli_SelecTripletsOriRel::cAppli_SelecTripletsOriRel
    (
        const std::vector<std::string> & aVArgs,
        const cSpecMMVII_Appli & aSpec,
        int aMode
        ) :
    cMMVII_Appli   (aVArgs,aSpec),
    mModeCompute   (aMode),
    mPhProj        (*this),
    mParamSel      (),
    mCandPool      (eTriCandPool::eAllOriTri),
    mNbMax         (mParamSel.mNbMax),
    mMinNb3        (5),
    mPropResFilter (0.85),
    mShow          (false)
{
}

cCollecSpecArg2007 & cAppli_SelecTripletsOriRel::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    if (mModeCompute==1)
    {
        anArgObl << Arg2007(mIm1,"name of image",{eTA2007::FileImage});
    }

    anArgObl  <<  mPhProj.DPOriRel().ArgDirInMand()
             <<  mPhProj.DPOrient().ArgDirInMand("Input orientation for calibration")  ;

    return anArgObl;
}

cCollecSpecArg2007 & cAppli_SelecTripletsOriRel::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    anArgOpt
        << cHeaderSectionArg("Input for tie points (at least one must  be provided)")
        << mPhProj.AddAllDpTieP(anArgOpt);

    anArgOpt << cHeaderSectionArg("Selection");
    if (mModeCompute==1)
        anArgOpt << AOpt2007(mIm2,"Im2","Process only the pair (Im1,Im2), print only (debug)");

    anArgOpt
        << AOpt2007(mCandPool,"CandPool","Pool of candidate third images",{eTA2007::HDV})
        << AOpt2007(mParamSel.mGreedy,"Greedy","If false, all valid candidates are selected (pass-through)",{eTA2007::HDV})
        << AOpt2007(mNbMax,"NbMax","Max number of triplets selected per pair",{eTA2007::HDV})
        << AOpt2007(mParamSel.mKappa,"Kappa","Cap of gain update in cells not covered by selected image",{eTA2007::HDV})
        << AOpt2007(mParamSel.mGainThr,"GainThr","Stop when best normalised gain is below",{eTA2007::HDV})
        << AOpt2007(mParamSel.mNbCells,"NbCells","Max number of cells per side of the grid",{eTA2007::HDV})
        << AOpt2007(mParamSel.mBsHLim,"BsHLim","Half-saturation of the B/H gain",{eTA2007::HDV})
        << AOpt2007(mParamSel.mAttenDens,"AttenDens","Saturation of the coverage term",{eTA2007::HDV})
        << AOpt2007(mMinNb3,"MinNb3","Min number of triple tie points for a candidate",{eTA2007::HDV})
        << AOpt2007(mPropResFilter,"PropResFilter","Proportion of tie points kept by the residual filter",{eTA2007::HDV})
        << AOpt2007(mShow,"Show","Print details",{eTA2007::HDV})
        ;
    return anArgOpt;
}

void cAppli_SelecTripletsOriRel::ReadOrientedPairs()
{
    tNameSet aSetIm;
    ReadFromFile(aSetIm,mPhProj.OriRel_NameAllImages(true));
    mVNamesIm = ToVect(aSetIm);

    for (const auto & aNameIm : mVNamesIm)
    {
        std::string aNameFile = mPhProj.OriRel_NameOriAllPairsOf1Image(aNameIm,true);
        if (! ExistFile(aNameFile))
            continue;
        std::vector<cCdtFinalPoseRel2Im> aVRes;
        ReadFromFile(aVRes,aNameFile);
        for (const auto & aRes : aVRes)
        {
            if (aRes.mVCdt.empty())
                continue;
            MMVII_INTERNAL_ASSERT_always(aRes.mIm1<aRes.mIm2,"ReadOrientedPairs : pair not ordered");
            std::vector<tPoseR> & aVPose = mMapPoses[tPairNames(aRes.mIm1,aRes.mIm2)];
            for (const auto & aCdt : aRes.mVCdt)
                aVPose.push_back(aCdt.mPose);
            mNeighbours[aRes.mIm1].insert(aRes.mIm2);
            mNeighbours[aRes.mIm2].insert(aRes.mIm1);
        }
    }
}

void cAppli_SelecTripletsOriRel::ReadCycleTriplets()
{
    for (const auto & aNameIm : mVNamesIm)
    {
        // the copy is created by the driver before it overwrites the original
        std::string aNameFile = mPhProj.OriRel_NameAllTripletsOf1Image(aNameIm,true,"Cyc."+GlobTaggedNameDefSerial());
        if (! ExistFile(aNameFile))
            aNameFile = mPhProj.OriRel_NameAllTripletsOf1Image(aNameIm,true);
        if (! ExistFile(aNameFile))
            continue;
        tSet3N aSet3;
        ReadFromFile(aSet3,aNameFile);
        std::vector<const cTripletName *> aV3;
        aSet3.PutInVect(aV3,false);
        for (const auto & a3 : aV3)
        {
            const auto & aN = a3->mNames;  // sorted
            mThirdsOfCycle[tPairNames(aN[0],aN[1])].insert(aN[2]);
            mThirdsOfCycle[tPairNames(aN[0],aN[2])].insert(aN[1]);
            mThirdsOfCycle[tPairNames(aN[1],aN[2])].insert(aN[0]);
        }
    }
}

bool cAppli_SelecTripletsOriRel::IsOriented(const std::string & aN1,const std::string & aN2) const
{
    return (aN1<aN2) ? (mMapPoses.count(tPairNames(aN1,aN2))!=0) : (mMapPoses.count(tPairNames(aN2,aN1))!=0);
}

tPoseR cAppli_SelecTripletsOriRel::PoseOfPair(const std::string & aN1,const std::string & aN2) const
{
    if (aN1<aN2)
        return mMapPoses.at(tPairNames(aN1,aN2)).at(0);
    return mMapPoses.at(tPairNames(aN2,aN1)).at(0).MapInverse();
}

const cSetHomogCpleIm & cAppli_SelecTripletsOriRel::FilteredHom(const std::string & aN1,const std::string & aN2)
{
    tPairNames aKey(aN1,aN2);
    auto anIt = mCacheHom.find(aKey);
    if (anIt!=mCacheHom.end())
        return anIt->second;

    cSetHomogCpleIm & aSetH = mCacheHom[aKey];
    int aNbInit;
    mPhProj.ReadHomolMultiSrce(aNbInit,aSetH,aN1,aN2);
    MMVII_INTERNAL_ASSERT_User(aNbInit!=0,eTyUEr::eUnClassedError,"No source for homologous points");
    FilterHomOnPropResidualAng
        (
            aSetH,
            mMapPoses.at(aKey),
            *mPhProj.InternalCalibFromImage(aN1),
            *mPhProj.InternalCalibFromImage(aN2),
            mPropResFilter
            );
    return aSetH;
}

std::vector<std::string> cAppli_SelecTripletsOriRel::Candidates(const std::string & aNA,const std::string & aNB) const
{
    std::vector<std::string> aRes;
    if (mCandPool==eTriCandPool::eAllOriTri)
    {
        auto anIt = mNeighbours.find(aNA);
        if (anIt!=mNeighbours.end())
            for (const auto & aNC : anIt->second)
                if ((aNC!=aNB) && IsOriented(aNB,aNC))
                    aRes.push_back(aNC);
    }
    else
    {
        auto anIt = mThirdsOfCycle.find(tPairNames(aNA,aNB));
        if (anIt!=mThirdsOfCycle.end())
            for (const auto & aNC : anIt->second)
                if (IsOriented(aNA,aNC) && IsOriented(aNB,aNC))
                    aRes.push_back(aNC);
    }
    return aRes;
}

std::vector<cTripletName> cAppli_SelecTripletsOriRel::SelecOfPair(const std::string & aNA,const std::string & aNB)
{
    MMVII_INTERNAL_ASSERT_always(aNA<aNB,"SelecOfPair : names not ordered");
    std::vector<cTripletName> aRes;

    cPerspCamIntrCalib * aCalA = mPhProj.InternalCalibFromImage(aNA);
    cPerspCamIntrCalib * aCalB = mPhProj.InternalCalibFromImage(aNB);
    tPoseR aPoseAB = NormalizedPoseOfPair(PoseOfPair(aNA,aNB));

    // ----------- pair points : grid and depth -----------
    std::vector<cPt3dr> aVDirA;
    std::vector<cPt3dr> aVDirB;
    std::vector<cPt2dr> aVPts2;
    for (const auto & aCple : FilteredHom(aNA,aNB).SetH())
    {
        aVDirA.push_back(aCalA->DirBundle(aCple.mP1));
        aVDirB.push_back(aCalB->DirBundle(aCple.mP2));
        aVPts2.push_back(GridCoordOfDir(aVDirA.back()));
    }
    tREAL8 aDepthH = (aVPts2.empty()) ? -1.0 : DepthOfPair(aPoseAB,aVDirA,aVDirB);
    if (mShow)
        StdOut() << "=== Pair " << aNA << " " << aNB << " Nb2=" << aVPts2.size() << " H=" << aDepthH << "\n";
    if (aDepthH<=0)
        return aRes;

    cSelecTripletsCoverage aSelec(mParamSel,aVPts2,aPoseAB.Tr(),aDepthH);

    // ----------- candidates -----------
    std::vector<std::string> aVNameOfId;  // name of C for each Id given to the selection
    for (const auto & aNC : Candidates(aNA,aNB))
    {
        std::vector<std::string> aV3 {aNA,aNB,aNC};
        std::sort(aV3.begin(),aV3.end());
        size_t aKA = std::find(aV3.begin(),aV3.end(),aNA) - aV3.begin();
        size_t aKB = std::find(aV3.begin(),aV3.end(),aNB) - aV3.begin();
        size_t aKC = std::find(aV3.begin(),aV3.end(),aNC) - aV3.begin();

        //  merge of the filtered pairs => triple tracks
        cMemoryInterfImportHom aMIIH;
        for (size_t aK1=0 ; aK1<3 ; aK1++)
            for (size_t aK2=aK1+1 ; aK2<3 ; aK2++)
                aMIIH.Add(FilteredHom(aV3[aK1],aV3[aK2]),aV3[aK1],aV3[aK2]);
        cComputeMergeMulTieP aMerge(aV3,&aMIIH);

        cPerspCamIntrCalib * aCalC = mPhProj.InternalCalibFromImage(aNC);
        std::vector<cPt3dr> aVDirA3;
        std::vector<cPt3dr> aVDirB3;
        std::vector<cPt3dr> aVDirC3;
        std::vector<cPt2dr> aVPts3;
        auto anIt = aMerge.Pts().find(tConfigIm{0,1,2});
        if (anIt!=aMerge.Pts().end())
        {
            const std::vector<cPt2dr> & aVPIm = anIt->second.mVPIm;
            for (size_t aKP=0 ; aKP+2<aVPIm.size() ; aKP+=3)
            {
                aVDirA3.push_back(aCalA->DirBundle(aVPIm.at(aKP+aKA)));
                aVDirB3.push_back(aCalB->DirBundle(aVPIm.at(aKP+aKB)));
                aVDirC3.push_back(aCalC->DirBundle(aVPIm.at(aKP+aKC)));
                aVPts3.push_back(GridCoordOfDir(aVDirA3.back()));
            }
        }

        tPoseR aPoseAC = PoseOfPair(aNA,aNC);
        cPt3dr aCenterC;
        bool aOk = (aVPts3.size()>=size_t(mMinNb3))
                   && CenterOfThirdImage(aCenterC,aPoseAB,aPoseAC,aVDirA3,aVDirB3,aVDirC3,size_t(mMinNb3));
        if (! aOk)
        {
            if (mShow)
                StdOut() << "  C=" << aNC << " Nb3=" << aVPts3.size() << " REJECTED\n";
            continue;
        }

        cCandTriCov aCand;
        aCand.mId     = int(aVNameOfId.size());
        aCand.mCenter = aCenterC;
        aCand.mVPts3  = aVPts3;
        aSelec.AddCandidate(aCand);
        aVNameOfId.push_back(aNC);

        if (mShow)
        {
            cSelecTripletsCoverage::cInfoCand anInfo = aSelec.InfoOfCandidate(aVNameOfId.size()-1);
            StdOut() << "  C=" << aNC << " Nb3=" << aVPts3.size()
                     << " Center=" << aCenterC
                     << " GamA=" << anInfo.mGammaA << " GamB=" << anInfo.mGammaB
                     << " Cells=" << anInfo.mNbCellCov << "/" << anInfo.mNbCells
                     << " MeanCov=" << anInfo.mMeanCover
                     << " Gain0=" << anInfo.mGain
                     << " AngCycle=" << AngleCycleRot(aPoseAB,aPoseAC,PoseOfPair(aNB,aNC)) << "\n";
        }

    }

    // ----------- greedy selection -----------
    for (const auto & aIdGain : aSelec.Select())
    {
        const std::string & aNC = aVNameOfId.at(aIdGain.first);
        aRes.push_back(cTripletName(aNA,aNB,aNC));
        if (mShow)
            StdOut() << "  SELECTED " << aNC << " Gain=" << aIdGain.second << "\n";
    }
    return aRes;
}

int cAppli_SelecTripletsOriRel::Exe()
{
    mPhProj.FinishInit();
    mParamSel.mNbMax = mNbMax;

    ReadOrientedPairs();
    if (mCandPool==eTriCandPool::eCycleSel)
        ReadCycleTriplets();

    if (mModeCompute==1)
    {
        bool aOnlyPair = IsInit(&mIm2);
        if (aOnlyPair)
        {
            MMVII_INTERNAL_ASSERT_User(mIm1<mIm2,eTyUEr::eUnClassedError,"Im2 must be > Im1 (pair stored in folder of Im1) : "+mIm1+" "+mIm2);
            MMVII_INTERNAL_ASSERT_User(IsOriented(mIm1,mIm2),eTyUEr::eUnClassedError,"Pair not oriented : "+mIm1+" "+mIm2);
        }

        tSet3N aSet3;
        size_t aNbPair = 0;
        for (const auto & aKeyPoses : mMapPoses)
        {
            const tPairNames & aPair = aKeyPoses.first;
            if (aPair.first!=mIm1)                          // pairs stored in folder of Im1 : (Im1,B), Im1<B
                continue;
            if (aOnlyPair && (aPair.second!=mIm2))
                continue;
            aNbPair++;
            for (const auto & a3 : SelecOfPair(aPair.first,aPair.second))
                aSet3.Add(a3);
        }
        StdOut() << "Image " << mIm1 << " : " << aNbPair << " pairs, " << aSet3.size() << " triplets selected\n";

        if (! aOnlyPair)
        {
            mPhProj.DPOriRel().SetDirOutInIfNotInit();
            SaveInFile(aSet3,mPhProj.OriRel_NameSelecTripletsOf1Image(mIm1,false));
        }
        return EXIT_SUCCESS;
    }

    MMVII_INTERNAL_ERROR("cAppli_SelecTripletsOriRel : mode not implemented");
    return EXIT_FAILURE;
}

/* ====================================================== */
/*               OriPoseSelecTripletsOf1Im                */
/* ====================================================== */

tMMVII_UnikPApli Alloc_SelecTripletsOf1Im(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
    return tMMVII_UnikPApli(new cAppli_SelecTripletsOriRel(aVArgs,aSpec,1));
}

cSpecMMVII_Appli  TheSpec_SelecTripletsOf1Im
    (
        "OriPoseSelecTripletsOf1Im",
        Alloc_SelecTripletsOf1Im,
        "Select triplets of the oriented pairs of one image by image coverage and B/H",
        {eApF::Ori},
        {eApDT::TieP,eApDT::Orient},
        {eApDT::Orient},
        __FILE__
        );

}; // MMVII