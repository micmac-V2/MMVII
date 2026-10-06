//#include "MMVII_PCSens.h"

// Test commit

#include "cAppliNewFranges.h"




namespace MMVII
{

namespace NS_FrangesDetect
{


    /* =================================================== */
    /*      Overiding of  cMMVII_Appli                     */
    /* =================================================== */


cAppliNewFrange::cAppliNewFrange(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
    cMMVII_Appli      (aVArgs,aSpec),
    mPhProj           (*this),
    mZoomRed          (2.0),
    mDerFactZ1        (0.5),
    mSigmaTensZ1      (10.0),
    mDoSimul          (0),
    mNbVisuGen        (0),
    mHasMask          (false),
    mImMask           (cPt2di(1,1)),
    mImZ1             (cPt2di(1,1)),
    mDImZ1            (nullptr),
    mImRed            (cPt2di(1,1)),
    mDImRed           (nullptr),
    mImRedBlur        (cPt2di(1,1)),
    mDImRedBlur       (nullptr),
    mImGrad           (cPt2di(1,1)),
    mImTens           (cPt2di(1,1)),
    /*
    mImTx             (cPt2di(1,1)),
    mDImTx            (nullptr),
    mImTy             (cPt2di(1,1)),
    mDImTy            (nullptr),*/
    mImMaxHor         (cPt2di(1,1)),
    mDImMaxHor        (nullptr),
    mImMaxVert        (cPt2di(1,1)),
    mDImMaxVert       (nullptr),
    mImMaxLocTD         (cPt2di(1,1)),
    mDImMaxLocTD        (nullptr),
    mImTeta           (1),
    mDImTeta          (nullptr),
    mImTgt            (1),
    mDImTgt           (nullptr),
    mYC               (-1),
    mImIntegr         (1),
    mDImIntegr        (nullptr),
    mRadiomBackGround        (-1e9),
    mImRadFrange          (1),
    mDImRadFrange         (nullptr)
{
}

cAppliNewFrange::~cAppliNewFrange()
{
}

cCollecSpecArg2007 & cAppliNewFrange::ArgObl(cCollecSpecArg2007 & anArgObl)
{
      return    anArgObl
            <<  Arg2007(mPatImage,"Name of input Image", {eTA2007::FileDirProj,{eTA2007::MPatFile,"0"}})
      ;
}

cCollecSpecArg2007 & cAppliNewFrange::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    anArgOpt
              << AOpt2007(mZoomRed,"Zoom","Zoom Red process init" ,{eTA2007::HDV})
              << AOpt2007(mDerFactZ1,"DericheF","Factor of deriche filter (for Zoom=1)" ,{eTA2007::HDV})
              << AOpt2007(mSigmaTensZ1,"SigmaTens","Sigma for avaragin tensor (Z=1)" ,{eTA2007::HDV})
              << AOpt2007(mDoSimul,"DoSimul","Make simulation/syntheic images : 1 parab, 2 circles" ,{eTA2007::HDV})
              << AOpt2007(mPatVisu,"PatVisu","Pattern for generating visualization" )
              << mPhProj.DPMask().ArgDirInOpt()

            //  << AOpt2007(mSigCurv,"SigCurv","Sima for smoothig curve",{eTA2007::HDV})
            //  << AOpt2007(mIntY,"IntY","Interval for Y",{eTA2007::HDV})
            //  << AOpt2007(mNbIter,"NbIt","Number of iter initial smoothing",{eTA2007::HDV})
            //  << AOpt2007(mDoVisu,"DoVisu","Generate Visualisation ?",{eTA2007::HDV})
            //  << AOpt2007(mWidhMinAll,"MinWidth","Minimal witdh, general case",{eTA2007::HDV})
            //  << AOpt2007(mWidhMinBorder,"BorderMinWidth","Minimal witdh for border",{eTA2007::HDV})
            //  << AOpt2007(mMinHeightBorderRight,"MinHeightBorderRight","Minimal Heitgh for right border",{eTA2007::HDV})
            //  << AOpt2007(mParamStack,"Stack","Stacking parm [Pat,Nb,Mode] ",{{eTA2007::ISizeV,"[3,3]"}})
    ;

     return anArgOpt;
}


std::vector<cOneHelpSampleCmp>  cAppliNewFrange::Samples() const
{
   return
   {
       {"MMVII ExtractFranges Retiga_000000105.tif Sigma=[10,2] DoVisu=1"},
       {"MMVII ExtractFranges Retiga_.*.tif"}
   };
}

int cAppliNewFrange::Exe()
{
    mPhProj.FinishInit();

    if (RunMultiSet(0,0)) // Case several images in //
    {
       return ResultMultiSet();
    }
    DoOneImage(UniqueStr(0)); // Case 1 image (may be recalled from multiple)
    return EXIT_SUCCESS;
}

/* =================================================== */
/*              Specific functions                     */
/* =================================================== */





void cAppliNewFrange::ComputeRadiomCste()
{
    // Estimation of background, suppose to be constant; estimate at the center
    // where there is less franges

    {
        tREAL8 mPropEstBack = 1/3.0;

        std::vector<tREAL8>   aVRad;
        int aNbY=3;

        for (int anX=0 ; anX<mSzRed.x() ; anX++)
        {
            for (int aDy=-aNbY ; aDy<=aNbY ; aDy++)
                aVRad.push_back(mDImRed->GetV(cPt2di(anX,mYC+aDy)));
        }
        mRadiomBackGround = NC_KthVal(aVRad,mPropEstBack);
    }

    // Estimation of radiometry of frange, it's variable and a function of Y
    mImRadFrange = tIm1D(mSzRed.y());
    mDImRadFrange = &(mImRadFrange.DIm());

    for (int anY=0 ; anY<mSzRed.y() ; anY++)
    {
        std::vector<tREAL8>   aVRad;
        for (int anX=0 ; anX<mSzRed.x() ; anX++)
        {
            aVRad.push_back(mDImRedBlur->GetV(cPt2di(anX,anY)));
        }
        tREAL8 aVal = IKthVal(aVRad,aVRad.size()-10);
        // aVal = NC_KthVal(aVRad,1/3.0);
        mDImRadFrange->SetV(anY,aVal);
    }

    ExpFilterOfStdDev(*mDImRadFrange,5,50.0);
}



void  cAppliNewFrange::DoOneImage(const std::string & aNameIm)
{
    mNameIm = aNameIm;
    mHasMask = mPhProj.ImageHasMask(aNameIm);
    mImZ1 = tIm::FromFile(mNameIm);
    mDImZ1 = & (mImZ1.DIm());
    mImRed = mImZ1.BiCubicDeZoom(mZoomRed);
    mDImRed = &(mImRed.DIm());
    mSzRed = mDImRed->Sz();

    if (mHasMask)
    {
        MMVII_INTERNAL_ASSERT_User_UndefE(mZoomRed==(int)mZoomRed,"Non int DeZoom with mask");
        mImMask = cIm2D<tU_INT1>::FromFile(mPhProj.NameMaskOfImage(aNameIm));
        mImMask = mImMask.Decimate((int)mZoomRed);

        StdOut() << "SZ MASK= " << mImMask.DIm().Sz() << "\n";
    }

    if (mDoSimul)
        MakeImSimul();

    mImRedBlur = mImRed.Dup();
    mDImRedBlur = & (mImRedBlur.DIm());
    ExpFilterOfStdDev(*mDImRedBlur,5,2.0,5.0);

    DoTensorProcessing();
    ComputeRadiomCste();

    MakeImageMaxLoc();

    DoVisu();
}


};  // NS_FrangesDetect


using  namespace NS_FrangesDetect;

// ============================= Old version ===============

tMMVII_UnikPApli Alloc_NewFrangeDetect(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec)
{
      return tMMVII_UnikPApli(new cAppliNewFrange(aVArgs,aSpec));
}


cSpecMMVII_Appli  TheSpecAppliFranges_2
(
     "ExtractFranges_2",
      Alloc_NewFrangeDetect,
      "New new ... Extraction of Franges, version N (image filter)",
      {eApF::ImProc},
      {eApDT::Image},
      {eApDT::Console},
      __FILE__
);



}; //  namespace MMVII
