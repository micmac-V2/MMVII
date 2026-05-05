#include "cCodedTargetRefine.h"

namespace MMVII
{


    /**************************************************************************/
    /*
    * cCdTDiscr methods
    */
    /**************************************************************************/

    cCdTDiscr::cCdTDiscr(const std::string& aName, const std::string& aImName):
        mName   (aName),
        mImName (aImName),
        mExtent     (cPt2di(1,1)),
        mIm         (cPt2di(1,1)),
        mDIm        (nullptr),
        mCrop       (cPt2di(1,1)),
        mDCrop      (nullptr),
        mCano       (cPt2di(1,1)),
        mDCano      (nullptr),
        mSimu       (cPt2di(1,1)),
        mMasq       (cPt2di(1,1)),
        mDMasq      (nullptr)
        {
            //
        }

    tU_INT1 cCdTDiscr::CornersOnIm()
    {
        tU_INT1 aNb = 0;
        for (const auto& aPt : mVImCorners)
        {
            (void) aPt;//if (mDIm->Inside(ToI(aPt))){++aNb;}
        }
        return aNb;
    }

    tREAL8 cCdTDiscr::RansacTFOnBits()
    {
        //-> filter b/w bit centers
        std::vector<cPt2dr> aVBBits = {}, aVWBits = {};
        for (const cPt2dr& aBit : mVBitCenters)
        {
            if (mDCano->GetV(ToI(aBit)) == 0)
            {
                aVBBits.push_back(mCdT2Im.Value(aBit));
            } else {
                aVWBits.push_back(mCdT2Im.Value(aBit));
            }
        }
        //RansacTF();
        return 0;
    }

    tIm cCdTDiscr::Crop()
    {
        tIm aCropIm(mExtent.Sz());
        tDIm* aCropDIm = &aCropIm.DIm();

        aCropDIm->CropIn(mExtent.P0ByRef(), *mDIm);

        return aCropIm;
    }


    /*tIm cCdTDiscr::Sample()
    {
        tIm aSample = tIm(mFrame.Sz());//-> the image to fill with interpolated values

    }*/

    /*
     * Export methods
     */

    void cCdTDiscr::SaveCrop(const std::string& aDir)
    {
        tIm aCropIm = Crop();
        SaveIm(&aCropIm.DIm(), aDir + "Crop-CdT" + mName + '-' + mImName);
    }

    void cCdTDiscr::SaveMasq(const std::string& aDir)
    {
        SaveIm(&mMasq.DIm(), aDir + "Masq-CdT" + mName + '-' + mImName);
    }

    void cCdTDiscr::SaveIm(tDIm* aDIm, std::string aPath)
    {
        aDIm->ToFile(aPath);
    }

    /*
     * SET/GET
     */

    void cCdTDiscr::SetCdT2Im(cAff2D_r aCdT2Im)
    {
        mCdT2Im = aCdT2Im;
    }

    void cCdTDiscr::SetExtent(cRect2 aExt)
    {
        mExtent = aExt;
    }

    void cCdTDiscr::SetMasq(tIm aMasq)
    {
        mMasq = aMasq;
    }

    cRect2 cCdTDiscr::Extent()
    {
        return mExtent;
    }

    cPt2dr cCdTDiscr::CdT2Im(cPt2dr aPt, bool inverse)
    {
        return (inverse ? mCdT2Im.Inverse(aPt) : mCdT2Im.Value(aPt));
    }

    cAff2D_r cCdTDiscr::GetCdT2Im()
    {
        return mCdT2Im;
    }

    std::vector<cPt2dr> cCdTDiscr::VCdT2Im(std::vector<cPt2dr> aVPts, bool inverse)
    {
        std::vector<cPt2dr> aRes {};
        for (const auto& aPt : aVPts) aRes.push_back(CdT2Im(aPt, inverse));
        return aRes;
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
                << Arg2007(mFSpecName,"Xml/Json name for bit encoding struct",{{eTA2007::XmlOfTopTag,cFullSpecifTarget::TheMainTag}})
                << mPhProj.DPOrient().ArgDirInMand("Cameras absolute orientations")
                << mPhProj.DPGndPt2D().ArgDirInMand("Coded targets image measurements")
            ;
    }

    cCollecSpecArg2007 & cAppli_CodedTargetRefine::ArgOpt(cCollecSpecArg2007 & anArgOpt)
    {
        return anArgOpt
               << mPhProj.DPGndPt3D().ArgDirInOpt("Descr","CdT 3D Descriptions")
               << AOpt2007(mShow,"Show","Show useful details", {eTA2007::HDV})
               << AOpt2007(mVisu,"Visu","Save visualisation of results", {eTA2007::HDV})
            ;
    }

    cAppli_CodedTargetRefine::cAppli_CodedTargetRefine(const std::vector<std::string>& aVArgs,
                                                           const cSpecMMVII_Appli& aSpec):
        cMMVII_Appli(aVArgs, aSpec),
        mPhProj(*this),
        mVDescr({}),
        mIm(cPt2di(1,1)),
        mDIm(nullptr)
    {
        //
    }

    int cAppli_CodedTargetRefine::Exe()
    {
        //----- [0] Load project primitives
        mPhProj.FinishInit();
        std::vector<std::string> aVIm = VectMainSet(0);
        mFSpec.reset(cFullSpecifTarget::CreateFromFile(mFSpecName));

        //-> if descriptions are provided fill mVDescr
        if(!mPhProj.DPGndPt3D().DirInIsNONE())
        {
            ReadFromFile(mVDescr, cCdTDescr::NameFile(mPhProj, true));
        }

        //----- [A] DoOneImage

        for (const auto& aIm : aVIm)
        {
            mCam        = mPhProj.ReadCamPC(aIm, true);
            mIm         = tIm::FromFile(aIm);
            mSetImMes   = mPhProj.LoadMeasureIm(aIm);
            std::vector<cCdTDiscr>          aVDiscr;
            std::vector<cSaveExtrEllipe>    aVEll;

            ReadFromFile(aVEll, cSaveExtrEllipe::NameFile(mPhProj, mSetImMes, true));

            for (const auto& aEll : aVEll)
            {
                cCdTDiscr aDis = cCdTDiscr(aEll.mNameCode, aIm);
                BuildDiscr(aDis, aEll.mAffIm2Ref.MapInverse());
                aDis.SaveMasq(mPhProj.DirVisuAppli());
                //SampleDiscr
                //RefineDiscr
            }

            /*
             * Try do the same with not auto-extracted measurements
            if (!mVDescr.empty())
            {
                for (const auto& aDes : mVDescr)
                {
                    if (!aSetImMes.NameHasMeasure(aDes.mName))//-> avoid duplicate cCdTDiscr
                    {
                        cAff2D_r aCdT2Im = Descr2Aff(aDes, aCam);
                        aVDiscr.push_back(cCdTDiscr(aDes.mName, aImName, aCdT2Im, aIm, mFSpec));
                    }
                }
            }
            */

            //for (auto& aDis : aVDiscr)
            //{
            //    auto aNb = aDis.NbCornersOnCam();//-> computes number of CdT corners visible on aIm

            //    if (mVisu && aNb > 3)
            //    {
            //        aDis.SaveCrop(mPhProj.DirVisuAppli());
            //        //aDis.SaveMasq(mPhProj.DirVisuAppli());
            //    }
            //}
        }

        return EXIT_SUCCESS;
    }

    void cAppli_CodedTargetRefine::BuildDiscr(cCdTDiscr& aDis, cAff2D_r aCdT2Im)
    {

        aDis.SetCdT2Im(aCdT2Im);//-> set CdT/image mapping

        tIm aCdTIm                      = mFSpec->OneImTarget(*mFSpec->EncodingFromName(aDis.mName));
        tDIm* aDCdTIm                   = &aCdTIm.DIm();
        std::vector<cPt2dr> aVCdTCorn   = Corners(ToR(aDCdTIm->P0()), ToR(aDCdTIm->P1()));//-> CdT corners
        std::vector<cPt2dr> aVImCorn    = aDis.VCdT2Im(aVCdTCorn, true);//-> CdT image corners

        aDis.SetExtent(BBox(aVImCorn));//-> bounding box of CdT image corners
        tIm aMasq = tIm(aDis.Extent().Sz());
        tDIm* aDMasq  = &aMasq.DIm();

        StdOut() << "set masq " << aDis.mName << std::endl;
        for (const auto& aPix : aDis.Extent())
        {
            if (aDCdTIm->Inside(ToI(aDis.CdT2Im(ToR(aPix)))))
            {
                aDMasq->SetVTruncIfInside(aPix, 125);
            }
        }
        aDis.SetMasq(aMasq);
        StdOut() << "masq ok" << aDis.mName << std::endl;

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
     * Other useful methods/classes
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

    cRansacSol::cRansacSol(cPt2dr aSol, tREAL8 aL1Score):
        mSol(aSol),
        mL1Score(aL1Score)
    {
    }

    /*!
     * @brief RansacTF
     * @param aVBPts    : classified points as black values
     * @param aVWPts    : classified points as white values
     * @param aIm1      : in image
     * @param aIm2      : out image
     * @param aIt       : nb iterations (200)
     * @param aRDist    : minimal grey level distance between b/w values (50)
     * @param aDMasq     : pts to forget when computing solution score (nullptr)
     * @return
     */

    cRansacSol RansacTF(std::vector<cPt2dr> aVBPts, std::vector<cPt2dr> aVWPts, tIm& aIm1, tIm& aIm2,
                        int aIt, int aRDist, tDIm* aDMasq, tU_INT1 aMasqV)
    {
        //----- set primitives and load data images
        tREAL8              a = 1, b = 0;
        cPt2dr              aBestSol(1,0);
        int                 aBestL1 = 100000;
        tDIm* aDIm1 = &aIm1.DIm();
        tDIm* aDIm2 = &aIm2.DIm();

        //----- iterates same process *it* times to find the best (a,b) solution
        for (int ix=0;ix<aIt;++ix)
        {
            int aL1Score = 0;

            //-> choose 2 random pts from aVPts with distant grey level
            cPt2dr aBPt = aVBPts[RandUnif_N(aVBPts.size())], aWPt = aVWPts[RandUnif_N(aVWPts.size())];
            if (abs(aDIm1->GetV(ToI(aBPt)) - aDIm1->GetV(ToI(aWPt))) < aRDist) {--ix; continue;}

            //->get image value and compute (a,b)
            tU_INT1 aB1 = aDIm1->GetV(ToI(aBPt)), aB2 = aDIm1->GetV(ToI(aWPt));
            tU_INT1 aW1 = aDIm2->GetV(ToI(aBPt)), aW2 = aDIm2->GetV(ToI(aWPt));
            /*
             * aB2 = a * aB1 + b
             * aW2 = a * aW1 + b
             */
            a = (aB2 - aW2) / (aB1 - aW1);
            b = aB2 - a * aB1;

            //-> compute score of the solution
            for (const auto& aPix : *aDIm2)
            {
                if (aDMasq)//-> reject pixel based on masq
                {
                    if (aDMasq->GetV(aPix) == aMasqV){continue;}
                }
                tU_INT1 aVal = a * aDIm1->GetV(aPix) + b;
                aL1Score += aVal - aDIm2->GetV(aPix);
            }

            //-> change best score and sol and if we find better solution
            if (aL1Score < aBestL1)
            {
                aBestL1     = aL1Score;
                aBestSol    = cPt2dr(a,b);
            }
        }

        return cRansacSol(aBestSol, aBestL1);
    }

    cAff2D_r Descr2Aff(const cCdTDescr& aDes, cSensorCamPC* aCam)
    {
        const tREAL8& aR = aDes.mRes;
        std::vector<cPt2di> aVCorners = {cPt2di(0,0), cPt2di(aR,0), cPt2di(aR,aR), cPt2di(0,aR)};
        std::vector<cPt2dr> aVIn = {}, aVOut = {};

        for (const auto& aPt : aVCorners)
        {
            aVIn.push_back(ToR(aPt));
            aVOut.push_back(aCam->Ground2Image(aDes.mCdT2Gnd.Value(cPt3dr(aPt.x(), aPt.y(), 0))));
        }

        tREAL8      aRes;
        cAff2D_r    aAff;
        aAff = aAff.StdGlobEstimate(aVIn, aVOut, &aRes, nullptr, cParamCtrlOpt::Default());

        return aAff;
    }

    std::vector<cPt2dr> Corners(const cPt2dr& aP0, const cPt2dr& aP1)
    {
        return {aP0, cPt2dr(aP1.x(), aP0.y()), aP1, cPt2dr(aP0.x(), aP1.y())};
    }

}
