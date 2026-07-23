#include "cCodedTargetDescribe.h"

namespace MMVII
{
    const std::vector<cPt2dr> SqCorners = {cPt2dr(-1,-1), cPt2dr(1,-1), cPt2dr(1,1), cPt2dr(-1,1)};

/******************************************************************************/
/*
* cSetOfAugCdt
*/
/******************************************************************************/

    cSetOfAugCdt::cSetOfAugCdt():
        mVCdt({})
    {
    }

    void cSetOfAugCdt::Add(cAugCdt aCdt)
    {
        mVCdt.push_back(aCdt);
    }

    bool cSetOfAugCdt::NameHasAug(std::string aName)
    {
        bool aRes = false;
        for (const cAugCdt& aCdt : mVCdt)
        {
            if (aCdt.mName == aName)
            {
                aRes = true;
                break;
            }
        }
        return aRes;
    }

    cAugCdt* cSetOfAugCdt::CdtOfName(std::string aName)
    {
        cAugCdt* aRes = nullptr;
        for (cAugCdt& aCdt : mVCdt)
        {
            if (aCdt.mName == aName)
            {
                aRes = &aCdt;
            }
        }
        return aRes;
    }

    std::vector<cAugCdt>& cSetOfAugCdt::Cdts()
    {
        sort(mVCdt.begin(), mVCdt.end());
        return mVCdt;
    }

    std::string cSetOfAugCdt::NameFile(const cPhotogrammetricProject& aPhProj, bool aInput)
    {
        return  (aInput ? aPhProj.DPGndPt3D().FullDirIn() : aPhProj.DPGndPt3D().FullDirOut())
               + "Aug-"
               +  aPhProj.DPOrient().DirIn()
               + "."+ cMMVII_Appli::CurrentAppli().TaggedNameDefSerial();
    }

    void cSetOfAugCdt::AddData(const cAuxAr2007& anAuxParam)
    {
        cAuxAr2007 anAux("SetAugCdt",anAuxParam);
        MMVII::AddData(cAuxAr2007("Augmentations", anAux), mVCdt);
    }

    void AddData(const cAuxAr2007& anAux, cSetOfAugCdt& aSet)
    {
        aSet.AddData(anAux);
    }

/******************************************************************************/
/*
 * cExtract
 */
/******************************************************************************/

    cExtract::cExtract(const cSensorCamPC *aCam, cSaveExtrEllipe aEll):
        mCam (aCam),
        mEll (aEll)
    {
    }

/******************************************************************************/
/*
 * cAugCdt
 */
/******************************************************************************/

    cAugCdt::cAugCdt(std::string aName, std::shared_ptr<cFullSpecifTarget> aFSpec):
        mName (aName),
        mOKAug(false),
        mOKInter (false),
        mCenter (aFSpec->Center()),
        mFSpec (aFSpec)
    {
    }

    cAugCdt::cAugCdt():
        mFSpec (nullptr)
    {
    }

    bool cAugCdt::operator<(const cAugCdt& aAug) const
    {
        return mName < aAug.mName;
    }

    void cAugCdt::AddExtract(cExtract aExt)
    {
        mVExtracts.push_back(aExt);
    }

    tU_INT1 cAugCdt::NbExtracts() const
    {
        return mVExtracts.size();
    }

    std::vector<cPt2di> cAugCdt::Corners() const
    {
        std::vector<cPt2di> aRes = {};
        const cPt2dr& aC = mCenter;
        for (const auto& aP : SqCorners) aRes.push_back(ToI(aC + cPt2dr(aC.x()*aP.x(), aC.y()*aP.y())));
        return aRes;
    }

    std::vector<cPt3dr> cAugCdt::GndCorners() const
    {
        std::vector<cPt3dr> aRes = {};
        for (const auto& aP : Corners())
        {
            aRes.push_back(mRef2Gnd.Value(cPt3dr(aP.x(), aP.y(), 0)));
        }
        return aRes;
    }

    std::vector<cPt2dr> cAugCdt::GImCorners(cSensorCamPC* aCam)
    {
        std::vector<cPt2dr> aRes = {};
        tAff2Dr aRef2Glob = Ref2GImEstim(aCam);
        for (const auto& aPt : Corners())
        {
            aRes.push_back(aRef2Glob.Value(ToR(aPt)));
        }
        return aRes;
    }

    tIm cAugCdt::RefIm() const
    {
        tIm aRes(cPt2di(1,1));
        if (mFSpec != nullptr)
        {
            aRes = mFSpec->OneImTarget(*mFSpec->EncodingFromName(mName));
        }
        return aRes;
    }

    void cAugCdt::SetFSpec(std::shared_ptr<cFullSpecifTarget> aFSpec) {mFSpec = aFSpec;}

    void cAugCdt::Spatialize(tREAL8 aGndInterTol)
    {
        //----- init reference/ground correspondences
        std::vector<cPt3dr> aVMarks = {}, aVGndPts = {};
        cStdStatRes aStatRes;//-> for intersection validation

        for (const auto& aP : Corners())
        {
            aVMarks.push_back(cPt3dr(aP.x(), aP.y(), 0));//-> corners points planar coordinates and z=0
            std::vector<tSeg3dr> aVBundles = {};
            //----- intersection computation from all images extractions of current CdT
            for (const auto& aExt : mVExtracts)
            {
                cPt2dr aImP = aExt.mEll.mAffIm2Ref.Inverse(ToR(aP));//-> image point
                aVBundles.push_back(aExt.mCam->Image2Bundle(aImP));//-> camera-to-ground bundle
            }
            cPt3dr aGndP = BundleInters(aVBundles);//-> spatial intersection of bundles
            aVGndPts.push_back(aGndP);//-> spatial coordinates of the mark
            //----- intersection validation
            for (const auto& aExt : mVExtracts)
            {//-> ground distance between back projected ground corner and theoretical coordinates
                cPt2dr aBackProjP = aExt.mCam->Ground2Image(aGndP);
                cPt2dr aImP = aExt.mEll.mAffIm2Ref.Inverse(ToR(aP));
                tREAL8 aPixD = Norm2(aImP - aBackProjP);
                aStatRes.Add(aExt.mCam->Gen_GroundSamplingDistance(aGndP) * aPixD);
            }
            m3DPrec = aStatRes.Avg();
            if (m3DPrec <= aGndInterTol) mOKInter = true;
        }
        //----- estimation of spatial similarity from correspondences
        tREAL8 aRes;
        mRef2Gnd = mRef2Gnd.StdGlobEstimate(aVMarks, aVGndPts, &aRes, nullptr, cParamCtrlOpt::Default());
        //----- similarity validation -> smthg with aRes ??
        if (mOKInter) mOKAug = true;
    }

    cCdTDiscr cAugCdt::Discretize(cSensorCamPC* aCam, bool& isIn) const
    {
        //----- discretization = computing map between CdT theoretical image & glob image
        std::vector<cPt2dr> aVOut = {}, aVIn = {};
        for (const auto& aP : Corners())
        {
            aVIn.push_back(ToR(aP));
            aVOut.push_back(aCam->Ground2Image(mRef2Gnd.Value(cPt3dr(aP.x(), aP.y(), 0))));
        }
        tREAL8      aRes;
        tAff2Dr    aAff;
        aAff = aAff.StdGlobEstimate(aVIn, aVOut, &aRes, nullptr, cParamCtrlOpt::Default());
        isIn = aCam->IsVisibleOnImFrame(aAff.Value(mCenter));
        return cCdTDiscr(mName, aCam->NameImage(), aAff, true);
    }

    tAff2Dr cAugCdt::Ref2GImEstim(cSensorCamPC* aCam) const
    {
        //----- Ref2GIm Cdt reference image to global image affinity
        std::vector<cPt2dr> aVOut = {}, aVIn = {};
        for (const auto& aP : Corners())
        {
            aVIn.push_back(ToR(aP));
            aVOut.push_back(aCam->Ground2Image(mRef2Gnd.Value(cPt3dr(aP.x(), aP.y(), 0))));
        }
        tREAL8 aRes;
        tAff2Dr aAff;
        aAff = aAff.StdGlobEstimate(aVIn, aVOut, &aRes, nullptr, cParamCtrlOpt::Default());
        return aAff;
    }


    void cAugCdt::AddData(const cAuxAr2007& anAux)
    {
        MMVII::AddData(cAuxAr2007("Name", anAux), mName);
        MMVII::AddData(cAuxAr2007("Ref2Gnd", anAux), mRef2Gnd);
        MMVII::AddData(cAuxAr2007("Center", anAux), mCenter);
    }

    void AddData(const cAuxAr2007 &anAux, cAugCdt &anEx)
    {
        anEx.AddData(anAux);
    }

    std::string const cAugCdt::Show() const
    {
        return "CdT: " + mName + "\t | "
               + "Mul: " + std::to_string(NbExtracts()) + "\t | "
               + "Inter: " + (mOKInter ? "yes" : "no") + "\t | "
               + "Saved: " + (mOKAug ? "yes" : "no") + "\t | "
               + "3DPrec: " + (m3DPrec > 1e-6 ? std::to_string(m3DPrec) : "***") + '\n';
    }


/******************************************************************************/
/*
 * cAppli_CodedTargetDescribe
 */
/******************************************************************************/

    class cAppli_CodedTargetDescribe : public cMMVII_Appli
    {
    public:
        cAppli_CodedTargetDescribe(const std::vector<std::string>& aVArgs,
                                   const cSpecMMVII_Appli& aSpec);
    private:
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;
        void VisuAugCdT(const cSetOfAugCdt aOKAugSet, cSensorCamPC* aCam, std::string StrV="J");
        std::string NameVisu(const std::string & aIm, const std::string & aPref, const std::string aPost);
        cPhotogrammetricProject             mPhProj;
        std::string                         mSpecImIn;
        bool                                mShow;
        std::vector<std::string>            mVisu;
        std::string                         mFSpecName;
        std::shared_ptr<cFullSpecifTarget>  mFSpec;
        tREAL8                              mGndInterTol;
    };

    cCollecSpecArg2007& cAppli_CodedTargetDescribe::ArgObl(cCollecSpecArg2007& anArgObl)
    {
        return anArgObl
               << Arg2007(mSpecImIn, "Pattern/file of images", {{eTA2007::MPatFile,"0"}, {eTA2007::FileDirProj}})
               << Arg2007(mFSpecName,"Xml/Json name for bit encoding struct",{{eTA2007::XmlOfTopTag,cFullSpecifTarget::TheMainTag}})
               << mPhProj.DPOrient().ArgDirInMand("Cameras absolute orientations")
               << mPhProj.DPGndPt2D().ArgDirInMand("Coded targets image measurements")
               << mPhProj.DPGndPt3D().ArgDirOutMand("Output for coded targets description")
            ;
    }

    cCollecSpecArg2007 & cAppli_CodedTargetDescribe::ArgOpt(cCollecSpecArg2007 & anArgOpt)
    {
        return anArgOpt
               << AOpt2007(mGndInterTol,"GndInterTol","Ground tolerance for spatial intersection of corners", {eTA2007::HDV})
               << AOpt2007(mShow,"Show","Show useful details", {eTA2007::HDV})
               << AOpt2007(mVisu,"Visu", "Generate visualisation images with results [StrV(J-peg, T-iff), PatImV]", {{eTA2007::HDV, "['','*']"}})
            ;
    }

    cAppli_CodedTargetDescribe::cAppli_CodedTargetDescribe(const std::vector<std::string>& aVArgs,
                                                           const cSpecMMVII_Appli& aSpec):
        cMMVII_Appli(aVArgs, aSpec),
        mPhProj(*this),
        mShow (false),
        mVisu ({"","*"}),
        mFSpec (nullptr),
        mGndInterTol (1e-2)
    {
    }

    void cAppli_CodedTargetDescribe::VisuAugCdT(cSetOfAugCdt aOKAugSet, cSensorCamPC *aCam, std::string StrV)
    {
        StdOut() << "generation of: " << aCam->NameImage() << " visualization image\n";
        cRGBImage aIm = cRGBImage::FromFile(aCam->NameImage());
        auto aDIm = &aIm;
        aDIm->ResetGray();
        for (const auto& aAug : aOKAugSet.Cdts())
        {
            if (!aAug.mOKInter) continue;
            std::vector<cPt2dr> aVImP = {};
            auto& aCol = (aAug.mOKInter ? cRGBImage::Red : cRGBImage::Cyan);
            for (const auto& aP : aAug.GndCorners())
            {
                auto aImP = aCam->Ground2Image(aP);
                if (aDIm->InsideBL(aImP))
                {
                    aDIm->DrawCross(aImP, cPt2dr(1,1), aCol, .5);
                }
            }
        }
        std::string aNameV = NameVisu(aCam->NameImage(), "Marks", "");
        StrV == "J" ? aIm.ToJpgFileDeZoom(aNameV, 1, {"QUALITY=90"}) : aIm.ToFile(aNameV);
    }

    std::string cAppli_CodedTargetDescribe::NameVisu(const std::string & aIm, const std::string & aPref, const std::string aPost)
    {
        std::string aRes = mPhProj.DirVisuAppli() +  aPref +"-" + LastPrefix(FileOfPath(aIm));
        if (aPost!="") aRes = aRes + "-"+aPost;
        return aRes + ".tif";
    }

    int cAppli_CodedTargetDescribe::Exe()
    {
        //----- mPhProj prerequisites
        mPhProj.FinishInit();
        std::vector<std::string> aVIm = VectMainSet(0);
        mFSpec = std::shared_ptr<cFullSpecifTarget>(cFullSpecifTarget::CreateFromFile(mFSpecName));

        //----- parse augmentation candidates
        cSetOfAugCdt aAugSet;
        std::vector<std::string> aVNOriIm = {};//-> not oriented images

        for (const std::string& aIm : aVIm)
        {
            const cSensorCamPC* aCam = mPhProj.ReadCamPC(aIm, true, true);
            if (aCam == nullptr)//-> false if current image has been oriented
                {
                    aVNOriIm.push_back(aIm);//-> measurements not to be used
                    continue;//-> next image
                }

            //----- ellipses extrinsics extracted from current image
            std::vector<cSaveExtrEllipe> aVEll;
            ReadFromFile(aVEll, cSaveExtrEllipe::NameFile(mPhProj, mPhProj.LoadMeasureIm(aIm), true));

            for (const cSaveExtrEllipe& aEll : aVEll)
            {
                if (aAugSet.NameHasAug(aEll.mNameCode))//-> corresponding AugCdT already exists
                {
                    cAugCdt* aCdt = aAugSet.CdtOfName(aEll.mNameCode);
                    aCdt->AddExtract(cExtract(aCam, aEll));
                } else if (!starts_with(aEll.mNameCode, MMVII_NONE))//-> not found & not NONE target (=undecoded)
                {
                    cAugCdt aCdt(aEll.mNameCode, mFSpec);//-> new AugCdt
                    aCdt.AddExtract(cExtract(aCam, aEll));
                    aAugSet.Add(aCdt);
                }
            }
        }

        StdOut() << "-----\n"
                 << "rejected images (no orientation): " << aVNOriIm << '\n'
                << "--> " << aVNOriIm.size() << '/' << aVIm.size() << '\n'
                << "-----\n";

        //------ CdT spatialization
        cSetOfAugCdt aOKAugSet;//-> OK augmented CdT

        for (cAugCdt& aCdT : aAugSet.Cdts())
        {
            if (aCdT.NbExtracts() <= 3) continue;//-> no intersection with < 4 bundles
            aCdT.Spatialize(mGndInterTol);//-> spatial properties computation
            if (aCdT.mOKAug) aOKAugSet.Add(aCdT);
        }

        //----- show results & generate visualisation
        if (mShow) for (const auto& aCdT : aAugSet.Cdts()) StdOut() << aCdT.Show();
        StdOut() << "--> network augmentation: " << aOKAugSet.Cdts().size()
                 << '/' << aAugSet.Cdts().size() << '\n';

        if (mVisu[0] != "")
        {
            tNameSelector aSel = AllocRegex(mVisu[1]);
            for (const std::string& aIm : aVIm)
            {
                auto aCam = mPhProj.ReadCamPC(aIm, true, true);
                if (aCam == nullptr) continue;
                if (aSel.Match(aIm)) VisuAugCdT(aOKAugSet, aCam, mVisu[0]);
            }
        }

        //----- save result
        SaveInFile(aOKAugSet, cSetOfAugCdt::NameFile(mPhProj, false));//-> export to Aug-?.xml

        return EXIT_SUCCESS;
    }

    //----- memory allocation
    tMMVII_UnikPApli Alloc_CodedTargetDescribe(const std::vector<std::string> & aVArgs,
                                                  const cSpecMMVII_Appli & aSpec)
    {
        return tMMVII_UnikPApli(new cAppli_CodedTargetDescribe(aVArgs, aSpec));
    }

    cSpecMMVII_Appli TheSpec_CodedTargetDescribe
        (
            "CodedTargetDescribe",
            Alloc_CodedTargetDescribe,
            "Computes coded target 3D spatial properties from poses & images measurements",
            //metadonnees
            {eApF::Ori,eApF::GCP},//features
            {eApDT::ObjCoordWorld, eApDT::ObjMesInstr},//inputs
            {eApDT::Console},//output
            __FILE__
        );
}
