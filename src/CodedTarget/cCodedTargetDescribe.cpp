#include "cCodedTargetDescribe.h"

namespace MMVII
{
    const std::vector<cPt2dr> SqCorners = {cPt2dr(-1,-1), cPt2dr(1,-1), cPt2dr(1,1), cPt2dr(-1,1)};

    /**************************************************************************/
    /*
     * cExtract
     */
    /**************************************************************************/

    cExtract::cExtract(const cSensorCamPC *aCam, cSaveExtrEllipe aEll):
        mCam (aCam),
        mEll (aEll)
    {
    }

    /**************************************************************************/
    /*
     * cAugCdT
     */
    /**************************************************************************/

    cAugCdT::cAugCdT(std::string aName, std::shared_ptr<cFullSpecifTarget> aFSpec):
        mName (aName),
        mOKAug(false),
        mOKInter (false),
        mCenter (aFSpec->Center()),
        mFSpec (aFSpec)
    {
    }

    cAugCdT::cAugCdT():
        mFSpec (nullptr)
    {
    }

    bool cAugCdT::operator<(const cAugCdT& aAug) const
    {
        return mName < aAug.mName;
    }

    void cAugCdT::AddExtract(cExtract aExt)
    {
        mVExtracts.push_back(aExt);
    }

    tU_INT1 cAugCdT::NbExtracts() const
    {
        return mVExtracts.size();
    }

    std::vector<cPt2dr> cAugCdT::Corners() const
    {
        std::vector<cPt2dr> aRes = {};
        const cPt2dr& aC = mCenter;
        for (const auto& aP : SqCorners) aRes.push_back(aC + cPt2dr(aC.x()*aP.x(), aC.y()*aP.y()));
        return aRes;
    }

    std::vector<cPt3dr> cAugCdT::GndCorners() const
    {
        std::vector<cPt3dr> aRes = {};
        for (const auto& aP : Corners())
        {
            aRes.push_back(mRef2Gnd.Value(cPt3dr(aP.x(), aP.y(), 0)));
        }
        return aRes;
    }

    void cAugCdT::Spatialize(tREAL8 aGndInterTol)
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
                cPt2dr aImP = aExt.mEll.mAffIm2Ref.Inverse(aP);//-> image point
                aVBundles.push_back(aExt.mCam->Image2Bundle(aImP));//-> camera-to-ground bundle
            }
            cPt3dr aGndP = BundleInters(aVBundles);//-> spatial intersection of bundles
            aVGndPts.push_back(aGndP);//-> spatial coordinates of the mark
            //----- intersection validation
            for (const auto& aExt : mVExtracts)
            {//-> ground distance between back projected ground corner and theoretical coordinates
                cPt2dr aBackProjP = aExt.mCam->Ground2Image(aGndP);
                cPt2dr aImP = aExt.mEll.mAffIm2Ref.Inverse(aP);
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

    cCdTDiscr cAugCdT::Discretize(cSensorCamPC* aCam, bool& isIn) const
    {
        //----- discretization = computing map between CdT theoretical image & glob image
        std::vector<cPt2dr> aVOut;
        for (const auto& aP : Corners())
        {
            aVOut.push_back(aCam->Ground2Image(mRef2Gnd.Value(cPt3dr(aP.x(), aP.y(), 0))));
        }
        tREAL8      aRes;
        tAff2Dr    aAff;
        aAff = aAff.StdGlobEstimate(Corners(), aVOut, &aRes, nullptr, cParamCtrlOpt::Default());
        isIn = aCam->IsVisibleOnImFrame(aAff.Value(mCenter));
        return cCdTDiscr(mName, aCam->NameImage(), aAff);
    }

    void cAugCdT::AddData(const cAuxAr2007& anAux)
    {
        MMVII::AddData(cAuxAr2007("Name", anAux), mName);
        MMVII::AddData(cAuxAr2007("Ref2Gnd", anAux), mRef2Gnd);
        MMVII::AddData(cAuxAr2007("Center", anAux), mCenter);
    }

    void AddData(const cAuxAr2007 &anAux, cAugCdT &anEx)
    {
        anEx.AddData(anAux);
    }

    std::string cAugCdT::NameFile(const cPhotogrammetricProject& aPhProj, bool Input)
    {
        return  (Input ? aPhProj.DPGndPt3D().FullDirIn() : aPhProj.DPGndPt3D().FullDirOut())
               + "Aug-"
               +  aPhProj.DPOrient().DirIn()
               + "."+ cMMVII_Appli::CurrentAppli().TaggedNameDefSerial();
    }

    std::string cAugCdT::Show() const
    {
        return "CdT: " + mName + "\t | "
               + "Mul: " + std::to_string(NbExtracts()) + "\t | "
               + "Inter: " + (mOKInter ? "yes" : "no") + "\t | "
               + "Saved: " + (mOKAug ? "yes" : "no") + "\t | "
               + "3DPrec: " + (m3DPrec > 1e-6 ? std::to_string(m3DPrec) : "***") + '\n';
    }


    /**************************************************************************/
    /*
     * cAppli_CodedTargetDescribe
     */
    /**************************************************************************/

    class cAppli_CodedTargetDescribe : public cMMVII_Appli
    {
    public:
        cAppli_CodedTargetDescribe(const std::vector<std::string>& aVArgs,
                                   const cSpecMMVII_Appli& aSpec);
    private:
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;
        void VisuAugCdT(const std::vector<cAugCdT>& aVAugCdT, cSensorCamPC* aCam, std::string StrV="J");
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

    void cAppli_CodedTargetDescribe::VisuAugCdT(const std::vector<cAugCdT>& aVAugCdT, cSensorCamPC *aCam, std::string StrV)
    {
        StdOut() << "generation of: " << aCam->NameImage() << " visualization image\n";
        cRGBImage aIm = cRGBImage::FromFile(aCam->NameImage());
        auto aDIm = &aIm;
        aDIm->ResetGray();
        for (const auto& aAug : aVAugCdT)
        {
            if (!aAug.mOKInter) continue;
            std::vector<cPt2dr> aVImP = {};
            auto& aCol = (aAug.mOKInter ? cRGBImage::Red : cRGBImage::Cyan);
            //tU_INT1 ix = 0;
            for (const auto& aP : aAug.GndCorners())
            {
                //++ix;
                auto aImP = aCam->Ground2Image(aP);
                if (aDIm->InsideBL(aImP))
                {
                    aDIm->DrawCross(aImP, cPt2dr(1,1), aCol, .5);
                    //aDIm->DrawFiducial(aImP, cPt2dr(1,1), -1, aCol, .5);
                    //aDIm->DrawString(std::to_string(ix), cRGBImage::Green, aImP, cPt2dr(1,0));
                }
                //aVImP.push_back(aImP);
            }
            //auto& aC = mFSpec->Center();
            //auto aImC = aCam->Ground2Image(aAug.mRef2Gnd.Value(cPt3dr(aC.x(), aC.y(), 0)));
            //if (aDIm->InsideBL(aImC)) aDIm->DrawString(aAug.mName, aCol, aImC, cPt2dr(0.5,0.05));//-> name
            //aDIm->DrawPolygon(aVImP, cRGBImage::White, .5);//-> border
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

        //----- single image processing: init aVAugCdT
        std::vector<cAugCdT> aVAugCdT = {};
        std::vector<std::string> aVNOriIm = {};//-> not oriented images

        for (const std::string& aIm : aVIm)
        {
            const cSensorCamPC* aCam = mPhProj.ReadCamPC(aIm, true, true);
            if (aCam == nullptr)//-> checks if image has been oriented
                {
                    aVNOriIm.push_back(aIm);
                    continue;//-> reject image measurements
                }

            std::vector<cSaveExtrEllipe> aVEll;//-> load saved ellipses extrinsics from target extraction
            ReadFromFile(aVEll, cSaveExtrEllipe::NameFile(mPhProj, mPhProj.LoadMeasureIm(aIm), true));

            for (const cSaveExtrEllipe& aEll : aVEll)
            {
                bool isOK = false;
                for (cAugCdT& aCdT : aVAugCdT)//-> look if AugCdT has already been created
                {
                    if (aCdT.mName == aEll.mNameCode)//-> corresponding AugCdT already exists
                    {
                        aCdT.AddExtract(cExtract(aCam, aEll));//-> add extraction from current image
                        isOK = true;
                        break;
                    }
                }

                if (!isOK && !starts_with(aEll.mNameCode, MMVII_NONE))//-> if not found
                {
                    cAugCdT aCdT(aEll.mNameCode, mFSpec);//-> new AugCdT
                    aCdT.AddExtract(cExtract(aCam, aEll));//-> add extraction from current image
                    aVAugCdT.push_back(aCdT);
                }
            }
        }

        StdOut() << "-----\n"
                 << "rejected images (no orientation): " << aVNOriIm << '\n'
                << "--> " << aVNOriIm.size() << '/' << aVIm.size() << '\n'
                << "-----\n";

        sort(aVAugCdT.begin(), aVAugCdT.end());//-> sort vector on target names

        //------ CdT spatialization
        std::vector<cAugCdT> aVOKAugCdT = {};//-> OK augmented CdT

        for (cAugCdT& aCdT : aVAugCdT)
        {
            if (aCdT.NbExtracts() <= 3) continue;//-> no intersection with < 4 bundles
            aCdT.Spatialize(mGndInterTol);//-> spatial properties computation
            if (aCdT.mOKAug) aVOKAugCdT.push_back(aCdT);
        }

        //----- show results & generate visualisation
        if (mShow) for (const auto& aCdT : aVAugCdT) StdOut() << aCdT.Show();
        StdOut() << "--> network augmentation: " << aVOKAugCdT.size()
                 << '/' << aVAugCdT.size() << '\n';

        if (mVisu[0] != "")
        {
            tNameSelector aSel = AllocRegex(mVisu[1]);
            for (const std::string& aIm : aVIm)
            {
                auto aCam = mPhProj.ReadCamPC(aIm, true, true);
                if (aCam == nullptr) continue;
                if (aSel.Match(aIm)) VisuAugCdT(aVOKAugCdT, aCam, mVisu[0]);
            }
        }


        //----- save result
        SaveInFile(aVOKAugCdT, cAugCdT::NameFile(mPhProj, false));//-> export to Aug-?.xml

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
