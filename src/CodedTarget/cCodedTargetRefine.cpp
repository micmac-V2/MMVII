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
        mCdT        (cPt2di(1,1)),
        mSamp       (cPt2di(1,1)),
        mMask       (cPt2di(1,1))
        {
            //
        }

    tU_INT1 cCdTDiscr::CornersOnIm()
    {
        tU_INT1 aNb = 0;
        for (const auto& aPt : mVImCorners){if (mDIm->Inside(ToI(aPt))){++aNb;}}
        return aNb;
    }

    void cCdTDiscr::Sample()
    {
        mSamp = tIm(mExtent.Sz());
        tDIm* aDSamp = &mSamp.DIm();

        for (const auto& aPix : *aDSamp)
        {
            if (mMask.DIm().GetV(aPix) == MaskOutV)
            {
                //-> to avoid aliasing
                cRessampleWeigth aRW = cRessampleWeigth::GaussBiCub(ToR(aPix + mExtent.P0()), mCdT2Im.MapInverse(), 2);
                const std::vector<cPt2di>  & aVPts = aRW.mVPts;
                if (!aVPts.empty())
                {
                    tU_INT1 aV = 0;
                    for (int aK=0; aK<int(aVPts.size()) ; aK++)
                    {
                        if (mCdT.DIm().Inside(aVPts[aK]))
                        {
                            double aW = aRW.mVWeight[aK];
                            aV += aW * mCdT.DIm().GetV(aVPts[aK]);
                        }
                    }
                    aDSamp->SetV(aPix, aV);
                }
            }
        }
    }

    /*
     * export methods
     */

    void cCdTDiscr::SaveCrop(const std::string& aDir)
    {
        SaveIm(&mCrop.DIm(), aDir + "Crop-CdT" + mName + '-' + mImName);
    }

    void cCdTDiscr::SaveMask(const std::string& aDir)
    {
        SaveIm(&mMask.DIm(), aDir + "Mask-CdT" + mName + '-' + mImName);
    }

    void cCdTDiscr::SaveSample(const std::string& aDir)
    {
        SaveIm(&mSamp.DIm(), aDir + "Samp-CdT" + mName + '-' + mImName);
    }

    void cCdTDiscr::SaveIm(tDIm* aDIm, std::string aPath)
    {
        aDIm->ToFile(aPath);
    }

    void cCdTDiscr::SaveTmp(tIm& aTmp, const std::string& aDir, const std::string& aPref)
    {
        SaveIm(&aTmp.DIm(), aDir + aPref + mName + '-' + mImName);
    }

    /*
     * transformation methods
     */

    cPt2dr cCdTDiscr::CdT2Im(cPt2dr aPt, bool inverse) {return (inverse ? mCdT2Im.Inverse(aPt) : mCdT2Im.Value(aPt));}

    std::vector<cPt2dr> cCdTDiscr::VCdT2Im(std::vector<cPt2dr> aVPts, bool inverse)
    {
        std::vector<cPt2dr> aRes {};
        for (const auto& aPt : aVPts) aRes.push_back(CdT2Im(aPt, inverse));
        return aRes;
    }

    /*
     * set/get
     */

    cRect2  cCdTDiscr::Extent(){return mExtent;}
    tIm&    cCdTDiscr::CdT(){return mCdT;}
    tIm&    cCdTDiscr::Crop(){return mCrop;}
    tIm&    cCdTDiscr::Samp(){return mSamp;}
    tIm&    cCdTDiscr::Mask(){return mMask;}

    void    cCdTDiscr::SetCdT2Im(cAff2D_r aCdT2Im){mCdT2Im = aCdT2Im;}
    void    cCdTDiscr::SetExtent(cRect2 aExt){mExtent = aExt;}
    void    cCdTDiscr::SetCdT(tIm& aCdT){mCdT = aCdT;}
    void    cCdTDiscr::SetMask(tIm& aMask){mMask = aMask;}
    void    cCdTDiscr::SetCrop(tIm& aCrop){mCrop = aCrop;};

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
        cMMVII_Appli    (aVArgs, aSpec),
        mPhProj         (*this),
        mVDescr         ({}),
        mIm             (cPt2di(1,1)),
        mDIm            (nullptr)
    {
        //
    }

    int cAppli_CodedTargetRefine::Exe()
    {
        //----- PhProj primitives
        mPhProj.FinishInit();
        std::vector<std::string> aVIm = VectMainSet(0);
        mFSpec.reset(cFullSpecifTarget::CreateFromFile(mFSpecName));

        //-> if CdT descriptions are provided -> fill mVDescr (useful for CdT that have not been extracted)
        if(!mPhProj.DPGndPt3D().DirInIsNONE())
        {
            ReadFromFile(mVDescr, cCdTDescr::NameFile(mPhProj, true));
        }

        //----- independent process for each image

        for (const auto& aIm : aVIm)
        {
            mCam        = mPhProj.ReadCamPC(aIm, true);
            mIm         = tIm::FromFile(aIm);
            mSetImMes   = mPhProj.LoadMeasureIm(aIm);
            std::vector<cCdTDiscr>          aVDiscr;
            std::vector<cSaveExtrEllipe>    aVEll;

            ReadFromFile(aVEll, cSaveExtrEllipe::NameFile(mPhProj, mSetImMes, true));

            //----- refinement for CdT that have been extracted
            for (const auto& aEll : aVEll)
            {
                cCdTDiscr aDis = cCdTDiscr(aEll.mNameCode, aIm);
                BuildDiscr(aDis, aEll.mAffIm2Ref.MapInverse());//-> set prerequisites to use cCdTDiscr
                aDis.Sample();//-> creates simulated CdT
                if (mShow)
                {
                    aDis.SaveSample(mPhProj.DirVisuAppli());
                }
                DiscrMapRefine(aDis);
            }

            //----- do something close with some extra steps for unextracted CdT
            /*
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
        }

        return EXIT_SUCCESS;
    }

    void cAppli_CodedTargetRefine::BuildDiscr(cCdTDiscr& aDis, cAff2D_r aCdT2Im)
    {
        aDis.SetCdT2Im(aCdT2Im);//-> set CdT <-> image mapping

        tIm aCdT = mFSpec->OneImTarget(*mFSpec->EncodingFromName(aDis.mName));
        aDis.SetCdT(aCdT);

        tDIm* aDCdT = &aDis.CdT().DIm();
        std::vector<cPt2dr> aVCdTCorn   = Corners(ToR(aDCdT->P0()), ToR(aDCdT->P1()));//-> CdT corners
        std::vector<cPt2dr> aVImCorn    = aDis.VCdT2Im(aVCdTCorn);//-> CdT image corners

        aDis.SetExtent(BBox(aVImCorn));//-> bounding box of CdT image corners

        //----- set CdT croped image

        tIm aCrop(aDis.Extent().Sz());
        aDis.SetCrop(aCrop);

        tDIm* aCropDIm = &aDis.Crop().DIm();
        aCropDIm->CropIn(aDis.Extent().P0ByRef(), mIm.DIm());

        //----- set CdT in/out mask
        tIm aMask = tIm(aDis.Extent().Sz());
        tDIm* aDMask  = &aMask.DIm();

        for (const auto& aPix : aDis.Extent())
        {
            tU_INT1 aVal = aDCdT->Inside(ToI(aDis.CdT2Im(ToR(aPix), true))) ? MaskOutV : MaskInV;
            aDMask->SetV(aPix - aDis.Extent().P0(), aVal);
        }

        aDis.SetMask(aMask);

        if (mVisu)
        {
            aDis.SaveMask(mPhProj.DirVisuAppli());
            aDis.SaveCrop(mPhProj.DirVisuAppli());
        }
    }

    void cAppli_CodedTargetRefine::DiscrMapRefine(cCdTDiscr& aDis)
    {
        //----- residual mask computation
        std::vector<cPt2dr> aVSampWC = {}, aVSampBC = {};//-> w/b bits centers in mSamp
        tDIm* aDCdT = &aDis.CdT().DIm();
        tDIm* aDSamp = &aDis.Samp().DIm();

        for (const auto& aC : mFSpec->BitsCenters())
        {
            auto aSampC = aDis.CdT2Im(aC) - ToR(aDis.Extent().P0());
            aDCdT->GetV(ToI(aC)) < 125 ? aVSampBC.push_back(aSampC) : aVSampWC.push_back(aSampC);
        }

        cRansacSol aATF = RansacATF(aVSampBC, aVSampWC, aDSamp, &aDis.Crop().DIm(), &aDis.Mask().DIm(), MaskInV);

        if (mShow)
        {
            StdOut() << aDis.mName << "RansacATF solution -> " << aATF.mSol << "\n";
        }

        if (mVisu)//-> save corrected sample CdT as CorSamp
        {
            tIm aCSamp = tIm(aDSamp->Sz());//-> corrected sample
            tDIm* aDCSamp = &aCSamp.DIm();
            tIm aRes = tIm(aDSamp->Sz());//-> aATF residual image
            tDIm* aDRes = &aRes.DIm();
            tDIm* aDCrop = &aDis.Crop().DIm();

            //----- set w/b pixels on bit centers
            for (const auto& aP : aVSampBC) {aDCSamp->SetV(ToI(aP), 255);}
            for (const auto& aP : aVSampWC) {aDCSamp->SetV(ToI(aP), 0);}

            //----- set corrected pixels from ransac sol
            for (const auto& aP : *aDSamp)
            {
                aDCSamp->SetVTrunc(aP, aATF.mSol.x()*aDSamp->GetV(aP) + aATF.mSol.x());
                aDRes->SetVTrunc(aP, abs(aATF.mSol.x()*aDSamp->GetV(aP) + aATF.mSol.x() - aDCrop->GetV(aP)));
            }

            aDis.SaveTmp(aCSamp, mPhProj.DirVisuAppli(), "CorSamp-");
            aDis.SaveTmp(aRes, mPhProj.DirVisuAppli(), "ResSamp-");
        }
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

    cRansacSol RansacATF(std::vector<cPt2dr> aVBPts, std::vector<cPt2dr> aVWPts, tDIm* aDIm1, tDIm* aDIm2,
                         tDIm* aDMask, tU_INT1 aMaskInV, int aIt, int aRDist)
    {
        //----- set primitives and load data images
        tREAL8              a = 1.0, b = 0.0;
        cPt2dr              aBestSol(1,0);
        int                 aBestL1 = 10000000;

        //----- iterates same process *it* times to find the best (a,b) solution
        for (int ix=0;ix<aIt;++ix)
        {
            int aL1Score = 0;

            //-> choose 2 w/b random pts with suffisant grey level distance
            cPt2dr aBPt = aVBPts[RandUnif_N(aVBPts.size())], aWPt = aVWPts[RandUnif_N(aVWPts.size())];
            if (abs(aDIm1->GetV(ToI(aBPt)) - aDIm1->GetV(ToI(aWPt))) < aRDist) {continue;}

            //->get image value and compute (a,b)
            tREAL8 aB1 = aDIm1->GetV(ToI(aBPt)), aB2 = aDIm2->GetV(ToI(aBPt));//-> must be tREAL8 to compute a & b as tREAL8
            tREAL8 aW1 = aDIm1->GetV(ToI(aWPt)), aW2 = aDIm2->GetV(ToI(aWPt));
            /*
             * aB2 = a * aB1 + b
             * aW2 = a * aW1 + b
             */
            a = (aB2 - aW2) / (aB1 - aW1);
            b = aB2 - a * aB1;

            //-> compute score of the solution
            for (const auto& aPix : *aDIm2)
            {
                if (aDMask)//-> reject pixel if they are located on the mask
                {
                    if (aDMask->GetV(aPix) == aMaskInV){continue;}
                }
                tU_INT1 aVal = a * aDIm1->GetV(aPix) + b;
                aL1Score += abs(aVal - aDIm2->GetV(aPix));
            }

            //-> change best score and sol and if we've found a better solution
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
