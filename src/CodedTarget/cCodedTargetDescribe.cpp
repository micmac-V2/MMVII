#include "cCodedTargetDescribe.h"

namespace MMVII
{
    /**************************************************************************/
    /*
     * cCdTDet methods
     */
    /**************************************************************************/

    cCdTDet::cCdTDet(const cSensorCamPC* aCam, cMesIm1Pt aMes, cAff2D_r aAff2D):
        mCam(aCam),
        mMes(aMes),
        mIm2Ref (aAff2D)
    {
        //
    }

    /**************************************************************************/
    /*
     * cCdTDes methods
     */
    /**************************************************************************/

    cCdTDes::cCdTDes(std::string aName, std::unique_ptr<cFullSpecifTarget>& aSpec):
        mName(aName),
        mEnc(aSpec->EncodingFromName(aName)),
        mRes(600),
        mVBitCenters2D(aSpec->BitsCenters()), //assumes that target is a square
        mVCdTCorners({cPt2dr(0,0), cPt2dr(mRes,0), cPt2dr(mRes,mRes), cPt2dr(0,mRes)})
    {
        for (const auto& aPt : mVCdTCorners){mVCdtCorners3D.push_back(cPt3dr(aPt.x(), aPt.y(), 0));}
    }

    void cCdTDes::AddDetect(const cSensorCamPC* aCam, cMesIm1Pt aMes, cAff2D_r aAff2D)
    {
        cCdTDet aDet(aCam, aMes, aAff2D);
        mVDetects.push_back(aDet);
    }

    void cCdTDes::InterGndCorners(bool& aShow)
    {
        for (const cPt2dr& aCorn : mVCdTCorners)
        {
            std::vector<tREAL8> aVRes = {};
            cPt3dr aInter = CdT2GndByInter(aCorn, &aVRes);
            mVGndCorners.push_back(aInter);
            if (aShow)
            {
                StdOut() << "3D BUNDLE INTER" << aCorn << " -> " << aInter << ":\n";
                for (decltype(aVRes.size()) ix=0; ix<aVRes.size(); ++ix)
                {
                    StdOut() << mVDetects[ix].mCam->NameImage() << " -> " << aVRes[ix] << '\n';
                }
            }
        }
    }

    void cCdTDes::Estimate3DSimil(std::vector<cPt3dr>& aVInPts, std::vector<cPt3dr>& aVOutPts, bool& aShow)
    {
        tREAL8 aRes;
        mSimil3D = mSimil3D.StdGlobEstimate(aVInPts, aVOutPts, &aRes, nullptr, cParamCtrlOpt::Default());
        if (aShow)
        {
            StdOut() << "3D SIMIL ESTIMATE -> " << aRes << mSimil3D.Tr() << '\n';
        }
    }

    void cCdTDes::EstimateSimil3DOnCorners(bool& aShow)
    {
        Estimate3DSimil(mVCdtCorners3D, mVGndCorners, aShow);
    }

    cPt3dr cCdTDes::CdT2GndByInter(const cPt2dr& aPt, std::vector<tREAL8>* aVRes)
    {
        std::vector<tSeg3dr> aVBundles;
        for (const auto& aDet : mVDetects)
        {
            cPt2dr aImPt = aDet.mIm2Ref.Inverse(aPt);
            aVBundles.push_back(aDet.mCam->Image2Bundle(aImPt));
        }
        cPt3dr aInter = BundleInters(aVBundles);

        if(aVRes)
        {
            for (const auto& aDet : mVDetects)
            {
                tREAL8 aRes = Norm2(aDet.mIm2Ref.Inverse(aPt)-aDet.mCam->Ground2Image(aInter));
                aVRes->push_back(aRes);
            }
        }

        return aInter;
    }

    cPt2dr cCdTDes::Gnd2CdT(cPt3dr& aPt, const cCdTDet& aDet)
    {
        cPt2dr aImPt = aDet.mCam->Ground2Image(aPt);
        return aDet.mIm2Ref.Value(aImPt);
    }

    void cCdTDes::AddData(const cAuxAr2007& anAux)
    {
        MMVII::AddData(cAuxAr2007("Name", anAux), mName);
        MMVII::AddData(cAuxAr2007("CdT2Gnd", anAux), mSimil3D);
    }

    void AddData(const cAuxAr2007 &anAux, cCdTDes &aCdTDes)
    {
        aCdTDes.AddData(anAux);
    }

    std::string cCdTDes::NameFile(const cPhotogrammetricProject & aPhProj, std::string aName, bool Input)
    {
        return  (Input ? aPhProj.DPGndPt3D().FullDirIn() : aPhProj.DPGndPt3D().FullDirOut())
               + "CdTDescript-"
               +  aName
               + "."+ cMMVII_Appli::CurrentAppli().TaggedNameDefSerial();
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
               << mPhProj.DPOrient().ArgDirInMand("Camera absolute orientation")
               << mPhProj.DPGndPt2D().ArgDirInMand("Coded targets image measurements")
               << mPhProj.DPGndPt3D().ArgDirInMand("Coded targets ground coordinates")
               << mPhProj.DPGndPt3D().ArgDirOutMand("Output for coded target description<< mPhProj.DPOrient().ArgDirInMand()")
            ;
    }

    cCollecSpecArg2007 & cAppli_CodedTargetDescribe::ArgOpt(cCollecSpecArg2007 & anArgOpt)
    {
        return anArgOpt
               << AOpt2007(mShow,"Show","show some useful details", {eTA2007::HDV})//hdv = has default value
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

    void cAppli_CodedTargetDescribe::AddCdTDes(std::string aName, std::unique_ptr<cFullSpecifTarget>& aSpec)
    {
        cCdTDes aDes(aName, aSpec);
        mVCdTDes.push_back(aDes);
    }

    int cAppli_CodedTargetDescribe::Exe()
    {
        //----- [0] Load project primitives
        mPhProj.FinishInit();
        std::vector<std::string> aVIm = VectMainSet(0);
        mFSpec.reset(cFullSpecifTarget::CreateFromFile(mFSpecName));

        //----- [1] Load image measures
        mPhProj.LoadGCP3D();

        for (const std::string& aIm : aVIm)
        {
            const cSensorCamPC* aCam = mPhProj.ReadCamPC(aIm, true);
            cSetMesPtOf1Im aSetImMes = mPhProj.LoadMeasureIm(aIm);

            std::vector<cSaveExtrEllipe> aVEll;
            ReadFromFile(aVEll, cSaveExtrEllipe::NameFile(mPhProj, aSetImMes, true));

            for (const cSaveExtrEllipe& aEll : aVEll)
            {
                bool isOK = false;
                for (cCdTDes& aDes : mVCdTDes)
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
                    AddCdTDes(aEll.mNameCode, mFSpec);
                    mVCdTDes.back().AddDetect(aCam, aSetImMes.MeasuresOfName(aEll.mNameCode), aEll.mAffIm2Ref);
                }
            }
        }

        //------ [2] Intersect corners/centers/bits

        for (cCdTDes aDes : mVCdTDes)
        {
            if (mShow) StdOut() << "CdT -> " << aDes.mName << ":\n";
            aDes.InterGndCorners(mShow);//-> computes mVGndCorners
            aDes.EstimateSimil3DOnCorners(mShow);//-> computes m3DSimil
            SaveInFile(aDes, cCdTDes::NameFile(mPhProj, aDes.mName, false));
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
            "generate target 3D description from poses & images measurements",
            //metadonnees
            {eApF::Ori,eApF::GCP},//features
            {eApDT::ObjCoordWorld, eApDT::ObjMesInstr},//inputs
            {eApDT::Console},//output
            __FILE__
        );


    }
