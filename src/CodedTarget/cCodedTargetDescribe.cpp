#include "cCodedTargetDescribe.h"

namespace MMVII
{

    cCollecSpecArg2007& cAppli_CodedTargetDescribe::ArgObl(cCollecSpecArg2007& anArgObl)
    {
        return anArgObl
               << Arg2007(mSpecImIn, "Pattern/file of images", {{eTA2007::MPatFile,"0"}, {eTA2007::FileDirProj}})
               << Arg2007(mFSpecName,"Xml/Json name for bit encoding struct",{{eTA2007::XmlOfTopTag,cFullSpecifTarget::TheMainTag}})
               << mPhProj.DPGndPt2D().ArgDirInMand()
               << mPhProj.DPGndPt3D().ArgDirInMand()
               << mPhProj.DPOrient().ArgDirInMand()
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

    void cAppli_CodedTargetDescribe::AddDesCdT(std::string aName, std::unique_ptr<cFullSpecifTarget>& aSpec)
    {
        cDesCdT aDes(aName, aSpec);
        mVDesCdT.push_back(aDes);
    }

    cDetCdT::cDetCdT(const cSensorCamPC* aCam, cMesIm1Pt aMes, cAff2D_r aAff2D):
        mCam(aCam),
        mMes(aMes),
        mAff2D (aAff2D)
    {
        //
    }

    cDesCdT::cDesCdT(std::string aName, std::unique_ptr<cFullSpecifTarget>& aSpec):
        mName(aName),
        mEnc(aSpec->EncodingFromName(aName)),
        mRes(600), //assumes that target is a square
        mVCorners2D({cPt2dr(0,0), cPt2dr(mRes,0), cPt2dr(mRes,mRes), cPt2dr(0,mRes)})
    {
    }

    void cDesCdT::AddDetect(const cSensorCamPC* aCam, cMesIm1Pt aMes, cAff2D_r aAff2D)
    {
        cDetCdT aDet(aCam, aMes, aAff2D);
        mVDetects.push_back(aDet);
    }

    void cDesCdT::InterCorners(bool& show)
    {
        for (const cPt2dr& aCorn : mVCorners2D)
        {
            std::vector<tREAL8> aVRes = {};
            CdT2Gnd(aCorn, &aVRes);
            if (show)
            {
                StdOut() << "CdT : " + mName << "-> (bundle inter " << aCorn << ") RES: \n";
                for (decltype(aVRes.size()) ix=0; ix<aVRes.size(); ++ix)
                {
                    StdOut() << mVDetects[ix].mCam->NameImage() << " -> " << aVRes[ix] << '\n';
                }
            }
        }
    }

    cPt3dr cDesCdT::CdT2Gnd(const cPt2dr& aPt, std::vector<tREAL8>* aVRes)
    {
        std::vector<tSeg3dr> aVBundles;
        for (const auto& aDet : mVDetects)
        {
            cPt2dr aImPt = aDet.mAff2D.Inverse(aPt);
            aVBundles.push_back(aDet.mCam->Image2Bundle(aImPt));
        }
        cPt3dr aInter = BundleInters(aVBundles);

        for (const auto& aDet : mVDetects)
        {
            tREAL8 aRes = Norm2(aPt-Gnd2CdT(aInter, aDet));
            aVRes->push_back(aRes);
        }

        return aInter;
    }

    cPt2dr cDesCdT::Gnd2CdT(cPt3dr& aPt, const cDetCdT& aDet)
    {
        cPt2dr aImPt = aDet.mCam->Ground2Image(aPt);
        return aDet.mAff2D.Value(aImPt);
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
                for (cDesCdT& aDes : mVDesCdT)
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
                    AddDesCdT(aEll.mNameCode, mFSpec);
                    mVDesCdT.back().AddDetect(aCam, aSetImMes.MeasuresOfName(aEll.mNameCode), aEll.mAffIm2Ref);
                }
            }
        }

        //------ [2] Intersect corners/centers/bits

        for (cDesCdT aDes : mVDesCdT)
        {
            if (mShow) StdOut() << "CdT -> " << aDes.mName << ":\n";
            aDes.InterGndCorners(mShow);//-> computes mVGndCorners
            aDes.Estimate3DSimilOnCorners(mShow);//-> computes m3DSimil
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
