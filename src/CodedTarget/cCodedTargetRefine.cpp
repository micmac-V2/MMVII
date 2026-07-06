#include "cCodedTargetRefine.h"
#include "MMVII_Interpolators.h"

namespace MMVII
{

const tU_INT1 MaskOutV = 255, MaskInV = 0;//-> Val(aPix) = MaskOutV i.e aPix is out of the mask area

    /**************************************************************************/
    /*
    * cCdTDiscr methods
    */
    /**************************************************************************/

    cCdTDiscr::cCdTDiscr(const std::string& aName, const std::string& aImName, MMVII::tAff2Dr aAff):
        mName       (aName),
        mImName     (aImName),
        mExtent     (cPt2di(1,1)),
        mIm         (cPt2di(1,1)),
        mDIm        (nullptr),
        mCrop       (cPt2di(1,1)),
        mRef        (cPt2di(1,1)),
        mSamp       (cPt2di(1,1)),
        mMask       (cPt2di(1,1)),
        mInlMask    (cPt2di(1,1)),
        mRef2Glob   (aAff)
        {
        }

    void cCdTDiscr::Sample()
    {
        mSamp = tIm(mExtent.Sz());
        tDIm* aDSamp = &mSamp.DIm();
        aDSamp->InitCste(0);

        for (const auto& aPix : cRect2(aDSamp->Dilate(0)))
        {
            if (mMask.DIm().GetV(aPix) == MaskOutV)
            {
                //-> to avoid aliasing
                cRessampleWeigth aRW = cRessampleWeigth::GaussBiCub(ToR(aPix + mExtent.P0()), mRef2Glob.MapInverse(), 2);
                const std::vector<cPt2di>& aVPts = aRW.mVPts;
                if (!aVPts.empty())
                {
                    tREAL8 aV = 0;
                    tREAL8 aSW = 0;
                    for (int aK=0; aK<int(aVPts.size()) ; aK++)
                    {
                        if (mRef.DIm().Inside(aVPts[aK]))
                        {
                            double aW = aRW.mVWeight[aK];
                            aSW += aW;
                            aV += aW * mRef.DIm().GetV(aVPts[aK]);
                        }
                    }
                    aDSamp->SetV(aPix, aV/aSW);
                }
            }
        }
    }

    void cCdTDiscr::RansacTFSm2Cr()//-> computes mTFSm2Cr as a transfert function from mSamp to mCrop
    {
        std::vector<cPt2dr> aVSampBC = SampBitC('B');
        std::vector<cPt2dr> aVSampWC = SampBitC('W');

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
        aDCSamp->InitCste(255);

        //----- set corrected pixels from ransac solution
        for (const auto& aP : *aDSamp)
        {
            if (mMask.DIm().GetV(aP)==MaskOutV)
            aDCSamp->SetVTrunc(aP, mTFSm2Cr.mSol.x()*aDSamp->GetV(aP) + mTFSm2Cr.mSol.x());
        }

        //----- set w/b pixels bit centers to check location
        for (const auto& aP : SampBitC('B')) {aDCSamp->SetV(ToI(aP), 255);}
        for (const auto& aP : SampBitC('W')) {aDCSamp->SetV(ToI(aP), 0);}

        tIm aRes = ResTFSm2Cr();//-> TF residual image

        //----- save images
        SaveTmp(aCSamp, aDir, "CorSamp");
        SaveTmp(aRes, aDir, "ResSamp");
        SaveTmp(mTFSm2Cr.mInlMask, aDir, "InlMask");
    }

    void cCdTDiscr::LS10ParamSm2Cr()
    {
        tDIm* aDSamp =      &mSamp.DIm();
        tDIm* aDCrop =      &mCrop.DIm();
        tDIm* aDInlMask =   &mInlMask.DIm();

        cLS10PSys aSys = cLS10PSys(aDSamp, aDCrop, aDInlMask);
        aSys.Build();
        mSm2Cr = aSys.Solve();
    }

    void cCdTDiscr::VisuLS10ParamSm2Cr(const std::string& aDir)
    {
        LS10ParamSm2Cr();

        tDIm& aDSamp = mSamp.DIm();
        tDIm& aDMask = mMask.DIm();

        tIm aIm(aDSamp.Sz());
        tDIm& aDIm = aIm.DIm();
        aDIm.InitCste(255);
        std::unique_ptr<cDiffInterpolator1D> aInterpol;
        aInterpol.reset(cDiffInterpolator1D::AllocFromNames({"Linear"}));

        std::vector<tREAL8> aSol = mSm2Cr.mVP;

        for (const auto& aP : aDSamp)
        {
            if (aDMask.GetV(aP)==MaskOutV)
            {
                if (aDSamp.InsideInterpolator(*aInterpol, ToR(aP)))
                {
                    auto [aVG,aG] = aDSamp.GetValueAndGradInterpol(*aInterpol,ToR(aP));
                    auto& aDPx = aG.x();
                    auto& aDPy = aG.y();
                    tREAL8 aVS = aDSamp.GetV(aP);

                    tREAL8 aV = aVS * (aSol[0]*aP.x() + aSol[1]*aP.y() + aSol[2]) + aSol[3]
                                - aDPx * (aSol[4]*aP.x() + aSol[5]*aP.y() + aSol[6])
                                - aDPy * (aSol[7]*aP.x() + aSol[8]*aP.y() + aSol[9]);

                    aDIm.SetVTrunc(aP, aV);
                }
            }
        }
        SaveTmp(aIm, aDir, "LSCor");
    }

    tIm cCdTDiscr::MaskInCB(bool ext)
    {
        tIm aInCBMask(mSamp.DIm().Sz());
        tDIm& aDInCBMask = aInCBMask.DIm();
        aDInCBMask.InitCste(MaskInV);

        for (const auto& aP : aDInCBMask)
        {
            tREAL8 aExt = ext ? mRefCB.mBCD - mRefCB.mLmA : 0;
            tU_INT1 aV = InsideCB(Ref2Samp(ToR(aP), true), aExt) ? MaskOutV : MaskInV;
            aDInCBMask.SetV(aP, aV);
        }

        return aInCBMask;
    }

    tIm cCdTDiscr::MaskInCt(int aD)
    {
        tIm aInCtMask(mSamp.DIm().Sz());
        tDIm& aDInCtMask = aInCtMask.DIm();
        aDInCtMask.InitCste(MaskInV);

        auto aR = cRect2(mRef.DIm().Dilate(-aD));//-> mask rect in ref image

        for (const auto& aP : aDInCtMask)

        {
            if (aR.InsideBL(mRef2Glob.Inverse(ToR(aP + mExtent.P0()))))//-> if inside dilated rectangle (ie outside mask)
            {
                aDInCtMask.SetV(aP, MaskOutV);
            }
        }

        return aInCtMask;
    }


    bool cCdTDiscr::InsideCB(cPt2dr aP, tREAL8 ext)
    {
        return Norm2(aP - mRefCB.mC) <= mRefCB.mLmA + ext;
    }

    /*
     * export methods
     */

    void cCdTDiscr::SaveCrop(const std::string& aDir)
    {
        SaveIm(&mCrop.DIm(), aDir, "Crop");
    }

    void cCdTDiscr::SaveMask(const std::string& aDir)
    {
        SaveIm(&mMask.DIm(), aDir, "Mask");
    }

    void cCdTDiscr::SaveSample(const std::string& aDir)
    {
        SaveIm(&mSamp.DIm(), aDir, "Samp");
    }

    void cCdTDiscr::SaveRef(const std::string& aDir)
    {
        SaveIm(&mRef.DIm(), aDir, "Ref");
    }

    void cCdTDiscr::SaveIm(tDIm* aDIm, std::string aDir, std::string aPref)
    {
        aDIm->ToFile(aDir + "CdT" + mName + '-' + aPref + '-' + mImName);
    }

    void cCdTDiscr::SaveTmp(tIm& aTmp, const std::string& aDir, const std::string& aPref)
    {
        SaveIm(&aTmp.DIm(), aDir, aPref);
    }

    /*
     * transformation methods
     */

    cPt2dr cCdTDiscr::Ref2Im(cPt2dr aPt, bool inv) {return (inv ? mRef2Glob.Inverse(aPt) : mRef2Glob.Value(aPt));}

    std::vector<cPt2dr> cCdTDiscr::VRef2Im(std::vector<cPt2dr> aVPts, bool inv)
    {
        std::vector<cPt2dr> aRes {};
        for (const auto& aPt : aVPts) aRes.push_back(Ref2Im(aPt, inv));
        return aRes;
    }

    cPt2dr cCdTDiscr::Ref2Samp(cPt2dr aP, bool inv)
    {
        return (inv ? Ref2Im(aP + ToR(mExtent.P0()), inv) : Ref2Im(aP) - ToR(mExtent.P0()));
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
    tIm&        cCdTDiscr::InlMask(){return mInlMask;}
    cRansacSol  cCdTDiscr::TFSm2Cr(){return mTFSm2Cr;}
    cLS10PSol   cCdTDiscr::LSSm2Cr(){return mSm2Cr;}

    std::vector<cPt2dr> cCdTDiscr::SampBitC(char aCol)
    {
        std::vector<cPt2dr> aRes {};

        if (! contains("WB", aCol)) return aRes;

        std::vector<cPt2dr>& aVCenters = (aCol == 'W') ? mVWCenters : mVBCenters;
        for (const auto& aC : aVCenters)
        {
            aRes.push_back(Ref2Im(aC) - ToR(Extent().P0()));
        }

        return aRes;
    }

    cPt2dr cCdTDiscr::SampC()
    {
        return Ref2Im(mRefCB.mC) - ToR(Extent().P0());
    }

    tIm cCdTDiscr::ResTFSm2Cr()
    {
        tDIm& aDSamp = mSamp.DIm();
        tDIm& aDCrop = mCrop.DIm();
        tDIm& aDMask = mMask.DIm();
        tIm aRes(aDSamp.Sz());
        tDIm& aDRes = aRes.DIm();
        aDRes.InitCste(255);
        for (const auto& aP : aDSamp)
        {
            if (aDMask.GetV(aP) == MaskInV) continue;
            aDRes.SetVTrunc(aP, abs(mTFSm2Cr.mSol.x()*aDSamp.GetV(aP) + mTFSm2Cr.mSol.x() - aDCrop.GetV(aP)));
        }
        return aRes;
    }

    void cCdTDiscr::SetRef2Im(tAff2Dr aRef2Im){mRef2Glob = aRef2Im;}
    void cCdTDiscr::SetExtent(cRect2 aExt){mExtent = aExt;}
    void cCdTDiscr::SetRef(tIm aRef){mRef = aRef;}
    void cCdTDiscr::SetMask(tIm aMask){mMask = aMask;}
    void cCdTDiscr::SetInlMask(tIm aInlMask){mInlMask = aInlMask;}
    void cCdTDiscr::SetCrop(tIm aCrop){mCrop = aCrop;};

    void cCdTDiscr::SetWBCenters(const std::vector<cPt2dr>& aV)
    {
        tDIm* aDRef = &mRef.DIm();
        for (const auto& aC : aV)
        {
            aDRef->GetV(ToI(aC)) < 125 ? mVBCenters.push_back(aC) : mVWCenters.push_back(aC);
        }
    }

    void cCdTDiscr::SetCB(std::unique_ptr<cFullSpecifTarget>& aFSpec)
    {
        cPt2dr aMBW = aFSpec->CornerlEl_BW();//-> middle black to white corner
        cPt2dr aLWB = aFSpec->CornerlEl_WB();//-> lower white to black corner
        cPt2dr aBC  = aFSpec->BitsCenters()[0];//-> random bit center
        cPt2dr aC   = cPt2dr(aLWB.x(), aMBW.y());//-> center
        mRefCB = cCBParams(aC, Norm2(aC - aMBW), Norm2(aC - aLWB), Norm2(aC - aBC));
    }

    /**************************************************************************/
    /*
     * cSetOfCdTDiscr
     */
    /**************************************************************************/

    cSetOfCdTDiscr::cSetOfCdTDiscr(std::string aImName):
        mImName (aImName)
    {
    }

    void cSetOfCdTDiscr::Add(cCdTDiscr aDis)
    {
        mVDiscr.push_back(aDis);
    }

    bool cSetOfCdTDiscr::HasCdTName(std::string aCdTName)
    {
        for (const auto& aDis : mVDiscr)
        {
            if (aDis.mName == aCdTName) return true;
        }
        return false;
    }

    std::vector<cCdTDiscr>cSetOfCdTDiscr::CdTDiscretizations() const {return mVDiscr;}

    /**************************************************************************/
    /*
     * cAppli_CodedTargetRefine
     */
    /**************************************************************************/


    cCollecSpecArg2007& cAppli_CodedTargetRefine::ArgObl(cCollecSpecArg2007& anArgObl)
    {
        return anArgObl
                << Arg2007(mSpecImIn, "Pattern/file of images", {{eTA2007::MPatFile,"0"},eTA2007::FileImage})
                << Arg2007(mFSpecName,"Xml/Json name for bit encoding struct",{{eTA2007::XmlOfTopTag,cFullSpecifTarget::TheMainTag}})
                << mPhProj.DPGndPt2D().ArgDirInMand("Coded targets image measurements")
            ;
    }

    cCollecSpecArg2007 & cAppli_CodedTargetRefine::ArgOpt(cCollecSpecArg2007 & anArgOpt)
    {
        return anArgOpt
               << mPhProj.DPGndPt2D().ArgDirOutOptWithDef("Refine", "OutMes", "CdT 2D refined measurements")
               << mPhProj.DPGndPt3D().ArgDirInOpt("NetAug","CdT network augmentation")
               << mPhProj.DPOrient().ArgDirInOpt("OriAug","Absolute orientation -> mandatory if using network augmentation")
               << AOpt2007(mShow,"Show","Show useful details", {eTA2007::HDV})
               << AOpt2007(mVisu,"Visu","Save visualisation of results", {eTA2007::HDV})
               << AOpt2007(mMaskDil,"MaskDil","Dilate Ref image to filter inliers", {eTA2007::HDV})
            ;
    }

    cAppli_CodedTargetRefine::cAppli_CodedTargetRefine(const std::vector<std::string>& aVArgs,
                                                           const cSpecMMVII_Appli& aSpec):
        cMMVII_Appli    (aVArgs, aSpec),
        mPhProj         (*this),
        mVAugCdT         ({}),
        mIm             (cPt2di(1,1)),
        mDIm            (nullptr),
        mL1Lim          (20),
        mMaskDil        (0)
    {
        //
    }

    void cAppli_CodedTargetRefine::Visu(cSetMesPtOf1Im& aSet)
    {
        StdOut() << "generation of: " << mGlobImN << " visualization image\n";
        cRGBImage aIm = cRGBImage::FromFile(mGlobImN);
        auto aDIm = &aIm;
        aDIm->ResetGray();
        for (const auto& aMes : aSet.Measures())
        {
            auto& aCol = cRGBImage::Red;
            aDIm->DrawCross(aMes.mPt, cPt2dr(1,1), aCol, .5);
            aDIm->DrawString(aMes.mNamePt, cRGBImage::White, aMes.mPt, cPt2dr(1,0));
        }
        std::string aNameV = NameVisu(mGlobImN, "RefMes", "");
        std::string aStrV = "T";
        aStrV == "J" ? aIm.ToJpgFileDeZoom(aNameV, 1, {"QUALITY=90"}) : aIm.ToFile(aNameV);
    }


    int cAppli_CodedTargetRefine::Exe()
    {
        //----- PhProj primitives
        mPhProj.FinishInit();
        std::vector<std::string> aVIm = VectMainSet(0);
        mFSpec.reset(cFullSpecifTarget::CreateFromFile(mFSpecName));
        //----- if supplied, load coded targets augmentation
        if (!mPhProj.DPGndPt3D().DirInIsNONE() && mPhProj.DPOrient().DirInIsInit())
        {
            ReadFromFile(mVAugCdT, cAugCdT::NameFile(mPhProj, true));
        }
        //----- single image process

        for (const auto& aIm : aVIm)
        {
            mGlobImN = aIm;
            mIm = tIm::FromFile(mGlobImN);
            cSetOfCdTDiscr aSetOfDiscr(mGlobImN);//-> collection of image CdT discretizations
            cSetMesPtOf1Im aSet(mGlobImN);//-> to save final image measurements

            if (mShow) StdOut() << "(Im):" << mGlobImN <<'\n';

            //----- load image measurements obtained from standard image processing
            mSetImMes   = mPhProj.LoadMeasureIm(mGlobImN);
            std::vector<cSaveExtrEllipe> aVEll;
            ReadFromFile(aVEll, cSaveExtrEllipe::NameFile(mPhProj, mSetImMes, true));

            for (const auto& aEll : aVEll)
            {
                if (starts_with(aEll.mNameCode, MMVII_NONE)) continue;//-> reject undecoded targets
                aSetOfDiscr.Add(cCdTDiscr(aEll.mNameCode, aIm, aEll.mAffIm2Ref.MapInverse()));
            }
            //----- complete with augmentations if provided
            if (!mVAugCdT.empty() && !mPhProj.DPOrient().DirInIsNONE())//mPhProj.DPOrient().DirInIsNONE()
            {
                    cSensorCamPC* aCam = mPhProj.ReadCamPCFromFolder(mPhProj.DPOrient().DirIn(), aIm, true, true);
                if (aCam != nullptr)
                {
                    for (const auto& aAug : mVAugCdT)
                    {
                        if (aSetOfDiscr.HasCdTName(aAug.mName)) continue;//-> if already provided by ellipse extraction
                        bool isIn = false;//-> is CdT center visible on image
                        cCdTDiscr aDis = aAug.Discretize(aCam, isIn);
                        StdOut() << "ADD: CdT" << aDis.mName << '\n';
                        if (isIn) aSetOfDiscr.Add(aDis);
                    }
                }
            }

            //----- refinement for CdT which have been extracted
            for (auto& aDis : aSetOfDiscr.CdTDiscretizations())
            {
                    //-> set prerequisites to use cCdTDiscr
                    bool isOk = true;
                    BuildDiscr(aDis, isOk);
                    if (isOk)
                    {
                        //-> creation of sampled CdT
                        aDis.Sample();
                        if (mShow)
                        {
                            aDis.SaveSample(mPhProj.DirVisuAppli());
                        }
                        //-> refinement of image measurement based on CdT sampling
                        DiscrMapRefine(aDis);

                        //-> save refined measurement
                        cPt2dr aC = aDis.SampC();
                        cPt2dr aCorC = ToR(aDis.Extent().P0()) + aC + aDis.LSSm2Cr().mAff.Value(aC);

                        if (mShow)
                        {
                            StdOut() << "-------------------\n"
                                     << "CdT n°" << aDis.mName << '\n'
                                     << "-----\n"
                                     << "LS Params: " << aDis.LSSm2Cr().mVP << '\n'
                                     <<" refined mes :" << aC << "->" << aCorC << '\n';
                        }

                        cMesIm1Pt aMes(aCorC, aDis.mName, 1);
                        aSet.AddMeasure(aMes);
                    }
                }
            if (mVisu) Visu(aSet);
            mPhProj.SaveMeasureIm(aSet);
            }
        return EXIT_SUCCESS;
    }

    void cAppli_CodedTargetRefine::BuildDiscr(cCdTDiscr& aDis, bool &isOk)
    {
        aDis.SetCB(mFSpec);

        tIm aRef = mFSpec->OneImTarget(*mFSpec->EncodingFromName(aDis.mName));
        aDis.SetRef(aRef);

        tDIm* aDRef = &aDis.Ref().DIm();
        std::vector<cPt2dr> aVRefCorn   = Corners(ToR(aDRef->P0()), ToR(aDRef->P1()));//-> Ref corners
        std::vector<cPt2dr> aVImCorn    = aDis.VRef2Im(aVRefCorn);//-> Ref image corners

        //----- fails if not all CdT corners are in the image
        for (const auto& aC : aVImCorn)
        {
            if (!mIm.DIm().Inside(ToI(aC)))
            {
                StdOut() << aC;
                if (mShow) StdOut() << "Out!" << '\n';
                isOk = false;
                return;
            }
        }

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
        aDMask->InitCste(MaskInV);//-> by default all extent is masked

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
            aDis.SaveRef(mPhProj.DirVisuAppli());
        }
    }

    void cAppli_CodedTargetRefine::DiscrMapRefine(cCdTDiscr& aDis)
    {
        //----- ransac computation of Samp to Crop Transfert Function to filter outliers
        //mVisu ? aDis.VisuRansacTFSm2Cr(mPhProj.DirVisuAppli()) : aDis.RansacTFSm2Cr();
        //if (mShow) StdOut() << "CdT n°" << aDis.mName << " -> ransac TFSm2Cr solution " << aDis.TFSm2Cr().mSol << '\n';

        //aDis.SetInlMask(aDis.TFSm2Cr().mInlMask);

        //----- samp CB insidness to filter outliers
        //tIm aInCBMask = aDis.MaskInCB(true);
        //aDis.SetInlMask(aInCBMask);

        tIm aInCtMask = aDis.MaskInCt(mMaskDil);
        aDis.SetInlMask(aInCtMask);

        //----- least square computation of Samp to Crop 10-params mapping
        mVisu ? aDis.VisuLS10ParamSm2Cr(mPhProj.DirVisuAppli()) : aDis.LS10ParamSm2Cr();

        //----- pseudo-significativity test on radio bias
        //if (aDis.LSSm2Cr().mVP[3] > 80)
        //{
        //    tIm aInCBExtMask = aDis.MaskInCB(true);
        //    aDis.SetInlMask(aInCBExtMask);
        //    mVisu ? aDis.VisuLS10ParamSm2Cr(mPhProj.DirVisuAppli()) : aDis.LS10ParamSm2Cr();
        //}

        if (mVisu)
        {
            aDis.SaveTmp(aDis.InlMask(), mPhProj.DirVisuAppli(), "InlMask");
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
            "CdT refinement",
            {eApF::CodedTarget,eApF::ImProc},
            {eApDT::Image,eApDT::Xml},
            {eApDT::Console},//output
            __FILE__
            );

    /**************************************************************************/
    /*
     * Other useful methods/classes
     */
    /**************************************************************************/

    cCBParams::cCBParams()
    {
    }

    cCBParams::cCBParams(cPt2dr aC, tREAL8 aMLA, tREAL8 aLmA, tREAL8 aBCD):
        mC      (aC),
        mLMA    (aMLA),
        mLmA    (aLmA),
        mBCD    (aBCD)
    {
    }

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
            aDInlMask.InitCste(aMaskInV);//-> by default all pixels are outliers
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
                    aDInlMask.SetV(aPix, MaskOutV);//-> add pixel to the inlier mask with out-of-mask value MaskOutV
                    aL1Score.Add(aD);//-> add residual
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

    cLS10PSys::cLS10PSys(tDIm* aDIm1, tDIm* aDIm2, tDIm* aDInlMask, tU_INT1 aMaskInV):
        mDIm1(aDIm1),
        mDIm2(aDIm2),
        mDInlMask (aDInlMask),
        mMaskInV (aMaskInV),
        mSys(10)//-> LS system initialisation to a 10 parameters system
    {
    }

    void cLS10PSys::Build()
    {
        std::unique_ptr<cDiffInterpolator1D> aInterpol;
        aInterpol.reset(cDiffInterpolator1D::AllocFromNames({"Linear"}));

        //----- filling system with inliers observations
        for (const auto& aP : *mDIm1)
        {
            if (mDInlMask->GetV(aP)==MaskOutV)//-> only points that are out-of-mask value in inliers mask
            {
                if (mDIm2->InsideInterpolator(*aInterpol, ToR(aP)))//-> only points that are interpolable
                {
                    auto [aV,aG] = mDIm2->GetValueAndGradInterpol(*aInterpol,ToR(aP));//-> compute gradient to get partial derivative
                    tREAL8 aV1 = (tREAL8)mDIm1->GetV(aP);
                    tREAL8 aV2 = (tREAL8)mDIm2->GetV(aP);

                    auto& aPDx = aG.x();//-> interpolation gives x/y first order partial derivative
                    auto& aPDy = aG.y();

                    cDenseVect<tREAL8> aVEqObs({aV1*aP.x(),
                                                aV1*aP.y(),
                                                aV1,
                                                1,
                                                -aPDx*aP.x(),
                                                -aPDx*aP.y(),
                                                -aPDx,
                                                -aPDy*aP.x(),
                                                -aPDy*aP.y(),
                                                -aPDy});

                    mSys.PublicAddObservation(1, aVEqObs, aV2);
                }
            }
        }
    }

    cLS10PSol cLS10PSys::Solve()
    {
        return cLS10PSol(mSys.PublicSolve());
    }

    cLS10PSol::cLS10PSol(cDenseVect<tREAL8> aVParams):
        mVP (aVParams.ToStdVect())
    {
        mAff = tAff2Dr(cPt2dr(mVP[6], mVP[9]),//-> Tr
                        cPt2dr(mVP[4], mVP[7]),//-> Vx
                        cPt2dr(mVP[5], mVP[8]));//-> Vy
    }

    cLS10PSol::cLS10PSol()
    {
    }


    std::vector<cPt2dr> Corners(const cPt2dr& aP0, const cPt2dr& aP1)
    {
        return {aP0, cPt2dr(aP1.x(), aP0.y()), aP1, cPt2dr(aP0.x(), aP1.y())};
    }

    std::string cAppli_CodedTargetRefine::NameVisu(const std::string & aIm, const std::string & aPref, const std::string aPost)
    {
        std::string aRes = mPhProj.DirVisuAppli() +  aPref +"-" + LastPrefix(FileOfPath(aIm));
        if (aPost!="") aRes = aRes + "-"+aPost;
        return aRes + ".tif";
    }

}
