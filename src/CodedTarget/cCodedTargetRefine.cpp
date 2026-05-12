#include "cCodedTargetRefine.h"

namespace MMVII
{


    /**************************************************************************/
    /*
    * cCdTDiscr methods
    */
    /**************************************************************************/

    cCdTDiscr::cCdTDiscr(const std::string& aName, const std::string& aImName):
        mName       (aName),
        mImName     (aImName),
        mExtent     (cPt2di(1,1)),
        mIm         (cPt2di(1,1)),
        mDIm        (nullptr),
        mCrop       (cPt2di(1,1)),
        mRef        (cPt2di(1,1)),
        mSamp       (cPt2di(1,1)),
        mMask       (cPt2di(1,1))
        {
            //
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
                cRessampleWeigth aRW = cRessampleWeigth::GaussBiCub(ToR(aPix + mExtent.P0()), mRef2Im.MapInverse(), 2);
                const std::vector<cPt2di>  & aVPts = aRW.mVPts;
                if (!aVPts.empty())
                {
                    tU_INT1 aV = 0;
                    for (int aK=0; aK<int(aVPts.size()) ; aK++)
                    {
                        if (mRef.DIm().Inside(aVPts[aK]))
                        {
                            double aW = aRW.mVWeight[aK];
                            aV += aW * mRef.DIm().GetV(aVPts[aK]);
                        }
                    }
                    aDSamp->SetV(aPix, aV);
                }
            }
        }
    }

    void cCdTDiscr::RansacTFSm2Cr()//-> computes mTFSm2Cr as a transfert function from mSamp to mCrop
    {
        std::vector<cPt2dr> aVSampBC = SampC('B');
        std::vector<cPt2dr> aVSampWC = SampC('W');

        //----- computes ransac solution on w/b bit centers
        mTFSm2Cr = RansacATF(aVSampBC, aVSampWC, &mSamp.DIm(), &mCrop.DIm(), &mMask.DIm(), 10, MaskInV);
    }

    void cCdTDiscr::VisuRansacTFSm2Cr(const std::string& aDir)
    {
        RansacTFSm2Cr();

        //----- save corrected sample and ransac resiudals as CorSamp/ResSamp
        tDIm* aDSamp = &mSamp.DIm();

        tIm aCSamp = tIm(aDSamp->Sz());//-> corrected sample image
        tDIm* aDCSamp = &aCSamp.DIm();

        //----- set w/b pixels bit centers to check location
        for (const auto& aP : SampC('B')) {aDCSamp->SetV(ToI(aP), 255);}
        for (const auto& aP : SampC('W')) {aDCSamp->SetV(ToI(aP), 0);}

        //----- set corrected pixels from ransac solution
        for (const auto& aP : *aDSamp)
        {
            aDCSamp->SetVTrunc(aP, mTFSm2Cr.mSol.x()*aDSamp->GetV(aP) + mTFSm2Cr.mSol.x());
        }

        tIm aRes = ResTFSm2Cr();//-> TF residual image

        //----- save images
        SaveTmp(aCSamp, aDir, "CorSamp-");
        SaveTmp(aRes, aDir, "ResSamp-");
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

    cPt2dr cCdTDiscr::Ref2Im(cPt2dr aPt, bool inverse) {return (inverse ? mRef2Im.Inverse(aPt) : mRef2Im.Value(aPt));}

    std::vector<cPt2dr> cCdTDiscr::VRef2Im(std::vector<cPt2dr> aVPts, bool inverse)
    {
        std::vector<cPt2dr> aRes {};
        for (const auto& aPt : aVPts) aRes.push_back(Ref2Im(aPt, inverse));
        return aRes;
    }

    std::string cCdTDiscr::Cont()//-> fetch object actual content
    {
        std::string aRes;
        if (mCrop.DIm().Sz() != cPt2di(1,1))    aRes.push_back('C');
        if (mRef.DIm().Sz()  != cPt2di(1,1))    aRes.push_back('R');
        if (mMask.DIm().Sz() != cPt2di(1,1))    aRes.push_back('M');
        if (mSamp.DIm().Sz() != cPt2di(1,1))    aRes.push_back('S');
        return aRes;
    }

    bool cCdTDiscr::ReqCt(std::string aRCt)
    {
        bool aRes = true;
        std::string aCt = Cont();
        for (const auto& aC : aRCt)
        {
            aRes = contains(aCt, aC) && aRes;
        }
        return aRes;
    }

    /*
     * set/get
     */

    cRect2      cCdTDiscr::Extent(){return mExtent;}
    tIm&        cCdTDiscr::Ref(){return mRef;}
    tIm&        cCdTDiscr::Crop(){return mCrop;}
    tIm&        cCdTDiscr::Samp(){return mSamp;}
    tIm&        cCdTDiscr::Mask(){return mMask;}
    cRansacSol  cCdTDiscr::TFSm2Cr(){return mTFSm2Cr;}

    std::vector<cPt2dr> cCdTDiscr::SampC(char aC)
    {
        std::vector<cPt2dr> aRes {};

        if (! contains("WB",aC)) return aRes;

        std::vector<cPt2dr>& aVCenters = (aC == 'W') ? mVWCenters : mVBCenters;
        for (const auto& aC : aVCenters)
        {
            aRes.push_back(Ref2Im(aC) - ToR(Extent().P0()));
        }

        return aRes;
    }

    tIm cCdTDiscr::ResTFSm2Cr()
    {
        tDIm& aDSamp = mSamp.DIm();
        tDIm& aDCrop = mCrop.DIm();
        tIm aRes(aDSamp.Sz());
        tDIm& aDRes = aRes.DIm();
        for (const auto& aP : aDSamp)
        {
            aDRes.SetVTrunc(aP, abs(mTFSm2Cr.mSol.x()*aDSamp.GetV(aP) + mTFSm2Cr.mSol.x() - aDCrop.GetV(aP)));
        }
        return aRes;
    }


    void    cCdTDiscr::SetRef2Im(cAff2D_r aRef2Im){mRef2Im = aRef2Im;}
    void    cCdTDiscr::SetExtent(cRect2 aExt){mExtent = aExt;}
    void    cCdTDiscr::SetRef(tIm& aRef){mRef = aRef;}
    void    cCdTDiscr::SetMask(tIm& aMask){mMask = aMask;}
    void    cCdTDiscr::SetCrop(tIm& aCrop){mCrop = aCrop;};

    void    cCdTDiscr::SetWBCenters(const std::vector<cPt2dr>& aV)
    {
        tDIm* aDRef = &mRef.DIm();
        for (const auto& aC : aV)
        {
            aDRef->GetV(ToI(aC)) < 125 ? mVBCenters.push_back(aC) : mVWCenters.push_back(aC);
        }
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
        cMMVII_Appli    (aVArgs, aSpec),
        mPhProj         (*this),
        mVDescr         ({}),
        mIm             (cPt2di(1,1)),
        mDIm            (nullptr),
        mL1Lim (20)
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
                if (!starts_with(aEll.mNameCode, MMVII_NONE))
                {
                    //-> set prerequisites to use cCdTDiscr
                    cCdTDiscr aDis = cCdTDiscr(aEll.mNameCode, aIm);
                    BuildDiscr(aDis, aEll.mAffIm2Ref.MapInverse());
                    //-> creation of sampled CdT
                    aDis.Sample();
                    if (mShow)
                    {
                        aDis.SaveSample(mPhProj.DirVisuAppli());
                    }
                    //-> refinement of image measurement based on CdT sampling
                    DiscrMapRefine(aDis);
                }
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

    void cAppli_CodedTargetRefine::BuildDiscr(cCdTDiscr& aDis, cAff2D_r aRef2Im)
    {
        aDis.SetRef2Im(aRef2Im);//-> set Ref to image mapping

        tIm aRef = mFSpec->OneImTarget(*mFSpec->EncodingFromName(aDis.mName));
        aDis.SetRef(aRef);

        tDIm* aDRef = &aDis.Ref().DIm();
        std::vector<cPt2dr> aVRefCorn   = Corners(ToR(aDRef->P0()), ToR(aDRef->P1()));//-> Ref corners
        std::vector<cPt2dr> aVImCorn    = aDis.VRef2Im(aVRefCorn);//-> Ref image corners

        aDis.SetExtent(BBox(aVImCorn));//-> bounding box of Ref image corners

        aDis.SetWBCenters(mFSpec->BitsCenters());//-> set w/b bit centers w.r.t Ref images coordinates

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
            tU_INT1 aVal = aDRef->Inside(ToI(aDis.Ref2Im(ToR(aPix), true))) ? MaskOutV : MaskInV;
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
        //----- ransac computation of Samp to Crop Transfert Function
        mVisu ? aDis.VisuRansacTFSm2Cr(mPhProj.DirVisuAppli()) : aDis.RansacTFSm2Cr();
        if (mShow) StdOut() << aDis.mName << " -> ransac Sm2Cr solution " << aDis.TFSm2Cr().mSol << '\n';

        //----- set outlier mask from TF computation
        tIm aRes = aDis.ResTFSm2Cr();
        tDIm& aDRes = aRes.DIm();
        tIm aOMask(aDRes.Sz());//-> outlier mask
        tDIm& aDOMask = aOMask.DIm();
        for (const auto& aP : aDRes)
        {
            tU_INT1 aV= aDRes.GetV(aP) > mL1Lim ? MaskInV : MaskOutV;
            aDOMask.SetV(aP, aV);
        }
        //aDis.SaveTmp(aOMask, mPhProj.DirVisuAppli(), "OMask-");


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

    cRansacSol::cRansacSol():
        mInlMask(cPt2di(1,1))
    {
    }

    cRansacSol::cRansacSol(cPt2dr aSol):
        mSol(aSol),
        mInlMask(cPt2di(1,1))
    {
    }

    cRansacSol::cRansacSol(cPt2dr aSol, cStdStatRes aL1Score, tIm aInlMask):
        mSol(aSol),
        mL1Score(aL1Score),
        mInlMask(aInlMask)
    {
    }

    cRansacSol RansacATF(std::vector<cPt2dr> aVBPts, std::vector<cPt2dr> aVWPts, tDIm* aDIm1, tDIm* aDIm2,
                         tDIm* aDMask, tU_INT1 aGVT, tU_INT1 aMaskInV, int aIt, int aRDist)
    {
        //----- set primitives and load data images
        tREAL8              a = 1.0, b = 0.0;
        cRansacSol          aBestSol(cPt2dr(1,0));
        int                 aBestInl = 0;
        tIm                 aInlMask(aDIm1->Sz());
        tDIm&               aDInlMask = aInlMask.DIm();

        //----- iterates same process *it* times to find the best (a,b) solution
        for (int ix=0;ix<aIt;++ix)
        {
            int aInliers = 0;
            cStdStatRes aL1Score;

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
                tREAL8  aD = abs((tREAL8)aVal - (tREAL8)aDIm2->GetV(aPix));
                if (aD < aGVT)
                {
                    ++aInliers;
                    aDInlMask.SetV(aPix, MaskInV);
                    aL1Score.Add(aD);
                }
            }

            //-> change best sol
            if (aInliers > aBestInl)
            {
                aBestInl = aInliers;
                aBestSol = cRansacSol(cPt2dr(a,b), aL1Score, aInlMask);
            }
        }

        return aBestSol;
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
