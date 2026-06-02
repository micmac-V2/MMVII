#include "cCodedTargetDescribe.h"

namespace MMVII
{
    /**************************************************************************/
    /*
     * cCdTDetec methods
     */
    /**************************************************************************/

    cCdTDetec::cCdTDetec(const cSensorCamPC* aCam, cMesIm1Pt aMes, cAff2D_r aAff2D):
        mCam(aCam),
        mMes(aMes),
        mIm2Ref (aAff2D)
    {
        //
    }

    /**************************************************************************/
    /*
     * cCdTDescr methods
     */
    /**************************************************************************/

    cCdTDescr::cCdTDescr(std::string aName, std::unique_ptr<cFullSpecifTarget>& aSpec):
        mName           (aName),
        mEnc            (aSpec->EncodingFromName(aName)),
        mVBitCenters2D  (aSpec->BitsCenters())
    {
        cPt2di aSup(0,0), aInf(ToI(aSpec->Center()) * 2);
        mVCorners = {aSup, cPt2di(aInf.x(), aSup.y()),
                     aInf, cPt2di(aSup.x(), aInf.y())};
        mRes = aInf.x();
    }

    cCdTDescr::cCdTDescr():
        m3DOK (false)
    {
        //
    }

    void cCdTDescr::AddDetect(const cSensorCamPC* aCam, cMesIm1Pt aMes, cAff2D_r aAff2D)
    {
        cCdTDetec aDet(aCam, aMes, aAff2D);
        mVDetec.push_back(aDet);
    }

    int cCdTDescr::NbDetec(){return mVDetec.size();}

    void cCdTDescr::InterGndCorners(bool& aShow)
    {
        for (const cPt2di& aCorn : mVCorners)
        {
            std::vector<tREAL8> aVRes = {};
            cPt3dr aInter = CdT2GndByInter(aCorn, &aVRes);
            mVGndCorners.push_back(aInter);
            if (aShow)
            {
                StdOut() << "3D BUNDLE INTER" << aCorn << " -> " << aInter << ":\n";
                for (decltype(aVRes.size()) ix=0; ix<aVRes.size(); ++ix)
                {
                    StdOut() << mVDetec[ix].mCam->NameImage() << " -> " << aVRes[ix] << '\n';
                }
            }
        }
    }

    void cCdTDescr::EstimateCdT2Gnd(std::vector<cPt3dr>& aVInPts, std::vector<cPt3dr>& aVOutPts, bool& aShow)
    {
        tREAL8 aRes;
        mCdT2Gnd = mCdT2Gnd.StdGlobEstimate(aVInPts, aVOutPts, &aRes, nullptr, cParamCtrlOpt::Default());
        if (aShow)
        {
            StdOut() << "3D SIMIL ESTIMATE -> " << aRes << mCdT2Gnd.Tr() << '\n';
        }
    }

    void cCdTDescr::EstimateCdT2GndOnCorners(bool& aShow)
    {
        std::vector<cPt3dr> aVCornersZ0 = {};
        for (const auto& aC : mVCorners)
        {
            aVCornersZ0.push_back(cPt3dr(aC.x(), aC.y(), 0));
        }
        EstimateCdT2Gnd(aVCornersZ0, mVGndCorners, aShow);
    }

    cPt3dr cCdTDescr::CdT2GndByInter(const cPt2di& aPt, std::vector<tREAL8>* aVRes)
    {
        std::vector<tSeg3dr> aVBundles;
        for (const auto& aDet : mVDetec)
        {
            cPt2dr aImPt = aDet.mIm2Ref.Inverse(ToR(aPt));
            aVBundles.push_back(aDet.mCam->Image2Bundle(aImPt));
        }
        cPt3dr aInter = BundleInters(aVBundles);

        if(aVRes)
        {
            for (const auto& aDet : mVDetec)
            {
                tREAL8 aRes = Norm2(aDet.mIm2Ref.Inverse(ToR(aPt))-aDet.mCam->Ground2Image(aInter));
                aVRes->push_back(aRes);
            }
        }

        return aInter;
    }

    cPt2dr cCdTDescr::Gnd2CdT(cPt3dr& aPt, const cCdTDetec& aDet)
    {
        cPt2dr aImPt = aDet.mCam->Ground2Image(aPt);
        return aDet.mIm2Ref.Value(aImPt);
    }

    cPt3dr cCdTDescr::CdT2GndBySimil(cPt2di aPt)
    {
        return mCdT2Gnd.Value(cPt3dr(aPt.x(), aPt.y(), 0));
    }

    /*
     * Serialization methods
     */

    void cCdTDescr::AddData(const cAuxAr2007& anAux)
    {
        MMVII::AddData(cAuxAr2007("Name", anAux), mName);
        MMVII::AddData(cAuxAr2007("CdT2Gnd", anAux), mCdT2Gnd);
        MMVII::AddData(cAuxAr2007("Res", anAux), mRes);
    }

    void AddData(const cAuxAr2007 &anAux, cCdTDescr &aCdTDescr)
    {
        aCdTDescr.AddData(anAux);
    }

    std::string cCdTDescr::NameFile(const cPhotogrammetricProject & aPhProj, bool Input)
    {
        return  (Input ? aPhProj.DPGndPt3D().FullDirIn() : aPhProj.DPGndPt3D().FullDirOut())
               + "CdTDescript-"
               +  aPhProj.DPOrient().DirIn()
               + "."+ cMMVII_Appli::CurrentAppli().TaggedNameDefSerial();
    }

    /* new stuff cAugCdT & cExtract */

    cExtract::cExtract(const cSensorCamPC*& aCam, cSaveExtrEllipe aEll):
        mCam (aCam),
        mEll (aEll)
    {
        //
    }

    cAugCdT::cAugCdT(std::string aName, std::unique_ptr<cFullSpecifTarget>& aFSpec):
        mName (aName),
        mFSpec (aFSpec)
    {
        //
    }

    void cAugCdT::AddExtract(cExtract aExt)
    {
        mVExtracts.push_back(aExt);
    }

    std::vector<cPt2dr> cAugCdT::Corners()//-> specific to IGN Indoor CdT pattern
    {
        std::vector<cPt2dr> aRes = {};
        cPt2dr aMBW = mFSpec->CornerlEl_BW();//-> middle black to white corner
        cPt2dr aC   = cPt2dr(aMBW.x(), aMBW.x());//-> target centre
        for (const auto& aP : SqCorners) {aRes.push_back(aC * aP);}
        return aRes;
    }

    void cAugCdT::Spatialize()
    {
        //----- init reference/ground correspondences
        std::vector<cPt3dr> aRefMarks = {}, aVGndPts = {};

        for (const auto& aP : Corners())//-> space intersection of reference marks from extractions
        {
            aRefMarks.push_back(cPt3dr(aP.x(), aP.y(), 0));//-> corners points planar coordinates and z=0
            std::vector<tSeg3dr> aVBundles = {};
            for (const auto& aExt : mVExtracts)
            {
                cPt2dr aImP = aExt.mEll.mAffIm2Ref.Inverse(aP);//-> image point
                aVBundles.push_back(aExt.mCam->Image2Bundle(aImP));//-> camera to ground bundle
            }
            cPt3dr aGndP = BundleInters(aVBundles);
            aVGndPts.push_back(aGndP);
        }
        //----- estimate spatial similarity
        tREAL8 aRes;
        mRef2Gnd = mRef2Gnd.StdGlobEstimate(aRefMarks, aVGndPts, &aRes, nullptr, cParamCtrlOpt::Default());
        StdOut() << aRes;
        mOKAug = true;
    }


    /**************************************************************************/
    /*
     * cAppli_CodedTargetDescribe methods
     */
    /**************************************************************************/

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
               << AOpt2007(mShow,"Show","Show useful details", {eTA2007::HDV})
            ;
    }

    cAppli_CodedTargetDescribe::cAppli_CodedTargetDescribe(const std::vector<std::string>& aVArgs,
                                                           const cSpecMMVII_Appli& aSpec):
        cMMVII_Appli(aVArgs, aSpec),
        mPhProj(*this),
        mFSpec (nullptr)
    {
        //
    }

    void cAppli_CodedTargetDescribe::AddDescr(std::string aName, std::unique_ptr<cFullSpecifTarget>& aSpec)
    {
        cCdTDescr aDes(aName, aSpec);
        mVDescr.push_back(aDes);
    }

    int cAppli_CodedTargetDescribe::Exe()
    {
        //----- [0] Load project primitives
        mPhProj.FinishInit();
        std::vector<std::string> aVIm = VectMainSet(0);
        mFSpec.reset(cFullSpecifTarget::CreateFromFile(mFSpecName));

        //----- [1] Single image processing: filling mVDescr
        for (const std::string& aIm : aVIm)
        {
            const cSensorCamPC* aCam = mPhProj.ReadCamPC(aIm, true, true);
            if (aCam == nullptr)//->checks if image has been oriented
                {
                    if (mShow) StdOut() << aIm << " not oriented -> rejected\n";
                    continue;
                }

            cSetMesPtOf1Im aSetImMes = mPhProj.LoadMeasureIm(aIm);
            std::vector<cSaveExtrEllipe> aVEll;
            ReadFromFile(aVEll, cSaveExtrEllipe::NameFile(mPhProj, aSetImMes, true));

            for (const cSaveExtrEllipe& aEll : aVEll)
            {
                bool isOK = false;
                for (cCdTDescr& aDes : mVDescr)
                {
                    if (aDes.mName == aEll.mNameCode)
                    {
                        aDes.AddDetect(aCam, aSetImMes.MeasuresOfName(aEll.mNameCode), aEll.mAffIm2Ref);
                        isOK = true;
                        break;
                    }
                }
                if (!isOK && !starts_with(aEll.mNameCode, MMVII_NONE))
                {
                    AddDescr(aEll.mNameCode, mFSpec);
                    mVDescr.back().AddDetect(aCam, aSetImMes.MeasuresOfName(aEll.mNameCode), aEll.mAffIm2Ref);
                }
            }
        }

        //------ [2] Space intersection of corners & 3D similitude estimation
        std::vector<cCdTDescr> aVNotOKDescr = {};//-> for details

        for (decltype(mVDescr.size()) ix = 0; ix < mVDescr.size(); ++ix)
        {
            cCdTDescr& aDescr = mVDescr[ix];
            if (mShow) StdOut() << "CdT -> " << mVDescr[ix].mName << ":\n";
            if (aDescr.NbDetec() <= 3)
            {
                if (mShow) StdOut() << " not enough detections ("
                             << aDescr.NbDetec() << ") -> spatial description will not be computed\n";
                aVNotOKDescr.push_back(aDescr);
                continue;
            }
            aDescr.InterGndCorners(mShow);//-> computes mVGndCorners
            aDescr.EstimateCdT2GndOnCorners(mShow);//-> computes m3DSimil
            mVOKDescr.push_back(aDescr);
        }

        //----- [3] Save description if similiraty has been computed
        SaveInFile(mVOKDescr, cCdTDescr::NameFile(mPhProj, false));//-> export to CdTDescript-?.xml

        if (mShow)
        {
            StdOut() << "unsaved CdT descriptions: ";
            for (const auto& aDes : aVNotOKDescr)
            {
                StdOut() << aDes.mName << "-";
            }
            StdOut() << " -> OK saved: " << mVOKDescr.size() << '/' << mVDescr.size() << "\n";
        }

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
            "Generates coded target 3D description from poses & images measurements",
            //metadonnees
            {eApF::Ori,eApF::GCP},//features
            {eApDT::ObjCoordWorld, eApDT::ObjMesInstr},//inputs
            {eApDT::Console},//output
            __FILE__
        );
}
