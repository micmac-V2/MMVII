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
    mImLabel          (cPt2di(1,1)),
    mDImLabel         (nullptr),
    mImRegion           (cPt2di(1,1)),
    mDImRegion          (nullptr),
    mImCumNormT        (1),
    mDImCumNormT       (nullptr),
    mImTeta            (1),
    mDImTeta           (nullptr),
    mImTgt             (1),
    mDImTgt            (nullptr),
    mYC                (-1),
    mImIntegr          (1),
    mDImIntegr         (nullptr),
    mImSmoothTgt       (cPt2di(1,1)),
    mDImSmoothTgt      (nullptr),
    mImProjInt         (1),
    mDImProjInt        (nullptr),
    mRadiomBackGround   (-1e9),
    mImRadFrange        (1),
    mDImRadFrange       (nullptr),
    mVirtExtVert        (nullptr)
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


// std::vector<int> aV1;
// std::vector<int> aV2;



cCollecSpecArg2007 & cAppliNewFrange::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{


    anArgOpt
              << AOpt2007(mZoomRed,"Zoom","Zoom Red process init" ,{eTA2007::HDV})
              << AOpt2007(mDerFactZ1,"DericheF","Factor of deriche filter (for Zoom=1)" ,{eTA2007::HDV})
              << AOpt2007(mSigmaTensZ1,"SigmaTens","Sigma for avaragin tensor (Z=1)" ,{eTA2007::HDV})
              << AOpt2007(mDoSimul,"DoSimul","Make simulation/syntheic images : 1 parab, 2 circles" ,{eTA2007::HDV})
              << AOpt2007(mPatVisu,"PatVisu","Pattern for generating visualization" )

              << AOpt2007(mSigTensors ,"SigTens","Sigmas on tensors")

            //  << AOpt2007(aV1,"TestV1","Test vector int in once")
            //  << AOpt2007(aV2,"TestV2","Test vector int in once",{eTA2007::CanRepeat})


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

   //  StdOut() << "VVV " << aV1 << " " << aV2 << "\n";

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

void  cAppliNewFrange::SetLabelOn(const cPt2di & aPix,const std::vector<eTypePix> & aVLabel)
{
   for (const auto aLab : aVLabel)
       SetLabelOn(aPix,aLab);
}

static constexpr int LabelRegBorder= 0;
static constexpr int LabelRegUnInit = -1;
static constexpr int LabelRegMin= -2;


void cAppliNewFrange::MakeImageRegion()
{
    mImRegion =  cIm2D<tINT2> (mSzRed,nullptr,eModeInitImage::eMIA_Null);
    mDImRegion = & mImRegion.DIm();

    mDImRegion->InitCste(LabelRegUnInit);
    mDImRegion->InitBorder(LabelRegBorder);

    // Make barrier with min region
    for (const auto& aPix: *mDImRegion)
    {
        if (LabelIsOn(aPix,eTypePix::eMinTens))
            mDImRegion->SetV(aPix,LabelRegMin);
    }

    // init the vec point to propagate wiht max line
    std::vector<cPt3di> aVPCur;
    for (const auto & aVert : mGrReg.AllVertices())
    {
        const auto & aCC = aVert->Attr();

        for (const auto& aPt: aCC.mVPts)
        {
            mDImRegion->SetV(aPt,aCC.mNum);
            aVPCur.push_back(TP3z(aPt,aCC.mNum));
        }
    }
    const auto & aVP4 = Alloc4Neighbourhood();

    // do the propagation
    while (!aVPCur.empty())
    {
        std::vector<cPt3di> aNewVP;
        for (const auto & aP3 : aVPCur)
        {
            cPt2di aCurP2 = Proj(aP3);
            for (const auto & aN4 : aVP4)
            {
                cPt2di aNewP = aCurP2 + aN4;
                tINT2 & aV = mDImRegion->GetReference_V(aNewP);
                if ((aV==LabelRegUnInit) || (aV==LabelRegMin))
                {
                    if (aV==LabelRegUnInit)
                        aNewVP.push_back(TP3z(aNewP,aP3.z()));
                     aV = aP3.z();
                }
            }

        }
        aVPCur =    aNewVP;
    }

    //  costruct the edge of the graph by identify pixels wih different label
    for (const auto & aPix1 : *mDImRegion)
    {
        for (const auto & aN4 : aVP4)
        {
            cPt2di aPix2 = aPix1 + aN4;
            if (mDImRegion->Inside(aPix2))
            {
                 int aNum1 = mDImRegion->GetV(aPix1);
                 int aNum2 = mDImRegion->GetV(aPix2);
                 if (      (aNum1!=aNum2) && (aNum1>=0) && (aNum2>=0)
                       && (!LabelIsOn(aPix1,eTypePix::eMinTens))
                       && (!LabelIsOn(aPix2,eTypePix::eMinTens))
                     )
                 {
                     tVertexReg & aS1 = mGrReg.VertexOfNum(aNum1);
                     tVertexReg & aS2 = mGrReg.VertexOfNum(aNum2);

                     tEdgeReg* anE = mGrReg.AddEdge(aS1,aS2,cAOrientReg(),cAOrientReg(),cASymReg(),true);
                     FakeUseIt(anE);
                 }
            }
        }

    }
}


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
        cStdStatRes  aStat;
        std::vector<tREAL8>   aVRad;
        for (int anX=0 ; anX<mSzRed.x() ; anX++)
        {
            tREAL8 aRad = mDImRedBlur->GetV(cPt2di(anX,anY));
            aStat.Add(aRad);
            aVRad.push_back(aRad);
        }
      //  tREAL8 aVal = IKthVal(aVRad,aVRad.size()-10);
        tREAL8 aVal = aStat.Avg() + aStat.DevStd();
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

    mImLabel =  cIm2D<tU_INT4> (mSzRed,nullptr,eModeInitImage::eMIA_Null);
    mDImLabel = &mImLabel.DIm();



    if (mHasMask)
    {
        MMVII_INTERNAL_ASSERT_User_UndefE(mZoomRed==(int)mZoomRed,"Non int DeZoom with mask");
        mImMask = cIm2D<tU_INT1>::FromFile(mPhProj.NameMaskOfImage(aNameIm));
        mImMask = mImMask.Decimate((int)mZoomRed);

        StdOut() << "SZ MASK= " << mImMask.DIm().Sz() << "\n";
    }

    if (mDoSimul)    DoTensorProcessing();

        MakeImSimul();

    mImRedBlur = mImRed.Dup();
    mDImRedBlur = & (mImRedBlur.DIm());
    ExpFilterOfStdDev(*mDImRedBlur,5,2.0,5.0);

    DoTensorProcessing();
    ExtractTensorLines();

    AnalyseCC();

    MakeImageRegion();
    //GraphCompute();

    DoVisu();


    //DeleteAllAndClear(mVCMax);
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
