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

    void cAppli_CodedTargetDescribe::AddDesCdT(const cOneEncoding* aEnc)
    {
        cDesCdT aDes(aEnc);
        mVDesCdT.push_back(aDes);
    }

    cDetCdT::cDetCdT(const cSensorCamPC* aCam, cMesIm1Pt aMes, cAff2D_r aAff2D):
        mCam(aCam),
        mMes(aMes),
        mAff2D (aAff2D)
    {
        //
    }

    cDesCdT::cDesCdT(const cOneEncoding* aEnc):
        mName(aEnc->Name()),
        mEnc(aEnc)
    {
        //
    }

    void cDesCdT::AddDetect(const cSensorCamPC* aCam, cMesIm1Pt aMes, cAff2D_r aAff2D)
    {
        cDetCdT aDet(aCam, aMes, aAff2D);
        mVDetects.push_back(aDet);
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
                        aDes.AddDetect(aCam, aSetImMes.MeasuresOfName(aDes.mName), aEll.mAffIm2Ref);
                        isOK = true;
                        break;
                    }
                }
                if (!isOK && !starts_with(aEll.mNameCode, MMVII_NONE))
                {
                    AddDesCdT(mFSpec->EncodingFromName(aEll.mNameCode));
                }
            }
        }

        for (auto aDes : mVDesCdT)
        {
            StdOut() << aDes.mName << " " << aDes.mVDetects.size() <<'\n';
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
