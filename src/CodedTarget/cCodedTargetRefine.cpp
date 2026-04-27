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
        mVisib      (0)
    {
        mIm = cIm2D<tU_INT1>::FromFile(mCam->NameImage());
    }

    int cCdTDiscr::visibility()
    {
        CamExtent();
        if (mVisib > 0)
        {
            GenSimul();
        }
        return mVisib;
    }

    void cCdTDiscr::CamExtent()
    {
        std::vector<cPt3dr> aVCorner = mDescr.get3DCornersOnSimil();
        if (mCam->IsVisible(aVCorner[0]) && mCam->IsVisible(aVCorner[2]))
        {
            cPt2dr aUl = mCam->Ground2Image(aVCorner[0]), aLr = mCam->Ground2Image(aVCorner[2]);
            mExtent = cPixBox<2>(ToI(aUl), ToI(aLr));
            mCrop   = cIm2D<tU_INT1>(mExtent.Sz());
            mDCrop->CropIn(mExtent.P0ByRef(), *mDIm);
            mVisib  = 1;
        }
    }

    void cCdTDiscr::GenSimul()
    {

    }

    cPt2dr cCdTDiscr::CdT2Cam(cPt2di aPix)
    {
        cPt3dr aPt = CdT2Gnd(aPix);
        return mCam->Ground2Image(aPt);
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

        //----- [1] DoOneImage

        for (const auto& aIm : aVIm)
        {
            cSensorCamPC* aCam = mPhProj.ReadCamPC(aIm, true);
            for (const auto& aDes : mVDescr)
            {
                cCdTDiscr aDis(aDes, aCam);
                int aVis = aDis.visibility();//-> computes lots of things
                (void) aVis;
                //if (aVis > 3)
                //{
                //    ...
                //}
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
}
