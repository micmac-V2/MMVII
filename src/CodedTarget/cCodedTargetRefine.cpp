#include "cCodedTargetRefine.h"

namespace MMVII
{


    /**************************************************************************/
    /*
         * cCdTDiscr methods
         */
    /**************************************************************************/

    cCdTDiscr::cCdTDiscr(cCdTDescr aDes, cSensorCamPC* aCam):
        mName       (aDes.mName),
        mDescr      (aDes),
        mCam        (aCam),
        mExtent     (cPt2di(1,1)),
        mIm         (cPt2di(1,1)),
        mDIm        (nullptr),
        mCrop       (cPt2di(1,1)),
        mDCrop      (nullptr),
        mCano       (cPt2di(1,1)),
        mDCano      (nullptr),
        mSimu       (cPt2di(1,1)),
        mVisib      (0),
        mVCamCorner({})
    {
        mIm = cIm2D<tU_INT1>::FromFile(mCam->NameImage());
        mDIm = &mIm.DIm();
    }

    void cCdTDiscr::CamExtent()
    {
        CamCorners();
        mExtent = BBox(mVCamCorner);
        if (mCam->IsVisibleOnImFrame(ToR(mExtent.P0ByRef())) && mCam->IsVisibleOnImFrame(ToR(mExtent.P1ByRef())))
        {
            mCrop   = cIm2D<tU_INT1>(mExtent.Sz());
            mDCrop  = &mCrop.DIm();
            mDCrop->CropIn(mExtent.P0ByRef(), *mDIm);
            mVisib  = 1;
        }
    }

    void cCdTDiscr::CamCorners()
    {
        if (mVCamCorner.empty())
        {

        }
    }

    bool cCdTDiscr::FullExtentOnCam()
    {
        std::vector<cPt2dr> aVCorn = {};
        std::vector<cPt3dr> aVGndCorn = mDescr.getGndCornersOnSimil();
        for (const auto& aCorn : aVGndCorn)
        {
            aVCorn.push_back(mCam->Ground2Image(aCorn));
        }
    }

    void cCdTDiscr::GenSimul()
    {
        //localAff2D
        //
    }

    void cCdTDiscr::SaveCrop(const std::string& aDir)
    {
        SaveIm(mCrop, aDir + "Crop-" + mName + mCam->NameImage());
    }

    void cCdTDiscr::SaveIm(cIm2D<tU_INT1>& aIm, std::string aPath)
    {
        cDataIm2D<tU_INT1>* aDIm = &aIm.DIm();
        aDIm->ToFile(aPath);
    }

    /**************************************************************************/
    /*
     * cAppli_CodedTargetRefine methods
     */
    /**************************************************************************/

    cCollecSpecArg2007& cAppli_CodedTargetRefine::ArgObl(cCollecSpecArg2007& anArgObl)
    {
        return anArgObl
               << Arg2007(mSpecImIn, "Pattern/file of images", {{eTA2007::MPatFile,"0"}, {eTA2007::FileDirProj}})
            //<< Arg2007(mFSpecName,"Xml/Json name for bit encoding struct",{{eTA2007::XmlOfTopTag,cFullSpecifTarget::TheMainTag}})
               << mPhProj.DPOrient().ArgDirInMand("Cameras absolute orientations")
            //<< mPhProj.DPGndPt2D().ArgDirInMand("Coded targets image measurements")
               << mPhProj.DPGndPt3D().ArgDirInMand("Coded targets descriptions")
            //   << mPhProj.DPGndPt3D().ArgDirOutMand("Output for coded targets description")
            ;
    }

    cCollecSpecArg2007 & cAppli_CodedTargetRefine::ArgOpt(cCollecSpecArg2007 & anArgOpt)
    {
        return anArgOpt
               << AOpt2007(mShow,"Show","Show useful details", {eTA2007::HDV})
               << AOpt2007(mVisu,"Visu","Save visualisation of results", {eTA2007::HDV})
            ;
    }

    cAppli_CodedTargetRefine::cAppli_CodedTargetRefine(const std::vector<std::string>& aVArgs,
                                                           const cSpecMMVII_Appli& aSpec):
        cMMVII_Appli(aVArgs, aSpec),
        mPhProj(*this),
        mVDescr({})
    {
        //
    }

    int cAppli_CodedTargetRefine::Exe()
    {
        //----- [0] Load project primitives
        mPhProj.FinishInit();
        std::vector<std::string> aVIm = VectMainSet(0);
        ReadFromFile(mVDescr, cCdTDescr::NameFile(mPhProj, true));

        //----- [A] DoOneImage

        for (const auto& aIm : aVIm)
        {
            cSensorCamPC* aCam = mPhProj.ReadCamPC(aIm, true);
            for (const auto& aDes : mVDescr)
            {
                //----- [1] Is CdT visible on Cam ?
                    //----- (a) Is CdT full extent on Cam ?
                cCdTDiscr aDis(aDes, aCam);
                int aVis = aDis.visibility();//-> computes lots of things
                if (aVis > 0 && mVisu)
                {
                    aDis.SaveCrop(mPhProj.DirVisuAppli());
                }
            }
        }

        return EXIT_SUCCESS;
    }

    //----- memory allocation
    tMMVII_UnikPApli Alloc_CodedTargetRefine(const std::vector<std::string> & aVArgs,
                                               const cSpecMMVII_Appli & aSpec)
    {
        return tMMVII_UnikPApli(new cAppli_CodedTargetRefine(aVArgs, aSpec));
    }

    cSpecMMVII_Appli TheSpec_CodedTargetRefine
        (
            "CodedTargetRefine",
            Alloc_CodedTargetRefine,
            "blablabla",
            //metadonnees
            {eApF::Ori,eApF::GCP},//features
            {eApDT::ObjCoordWorld, eApDT::ObjMesInstr},//inputs
            {eApDT::Console},//output
            __FILE__
            );

    /**************************************************************************/
    /*
     * Other useful methods
     */
    /**************************************************************************/

    cPixBox<2> BBox(std::vector<cPt2dr> aVPts, int aMin, int aMax)
    {
        cPt2dr aSup(aMax,aMax), aInf(aMin,aMin);
        for (const auto& aPt : aVPts)
        {
            if (aPt.x() <= aSup.x()){aSup.x() = aPt.x();}
            if (aPt.y() <= aSup.y()){aSup.y() = aPt.y();}
            if (aPt.x() >= aInf.x()){aInf.x() = aPt.x();}
            if (aPt.y() >= aInf.y()){aInf.y() = aPt.y();}
        }
        return cPixBox<2>(ToI(aSup), ToI(aInf));
    }
}
