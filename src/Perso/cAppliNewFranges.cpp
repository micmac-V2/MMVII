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
    mDoSimul          (false),
    mDoVisu           (1),
    mImZ1             (cPt2di(1,1)),
    mDImZ1            (nullptr),
    mImRed            (cPt2di(1,1)),
    mDImRed           (nullptr),
    mImRedBlur        (cPt2di(1,1)),
    mDImRedBlur       (nullptr),
    mImTx             (cPt2di(1,1)),
    mDImTx            (nullptr),
    mImTy             (cPt2di(1,1)),
    mDImTy            (nullptr),
    mImMaxHor         (cPt2di(1,1)),
    mDImMaxHor        (nullptr),
    mImMaxVert        (cPt2di(1,1)),
    mDImMaxVert       (nullptr),
    mImMaxLoc         (cPt2di(1,1)),
    mDImMaxLoc        (nullptr),
    mImTeta           (1),
    mDImTeta          (nullptr),
    mImTgt            (1),
    mDImTgt           (nullptr),
    mYC               (-1),
    mImIntegr         (1),
    mDImIntegr        (nullptr),
    mImVisu           (cPt2di(1,1)),
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
              << AOpt2007(mDoSimul,"DoSimul","Make simulation/syntheic images" ,{eTA2007::HDV})
              << AOpt2007(mDoVisu,"DoVisu","Level of visualisation generated ?" ,{eTA2007::HDV})

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


// The simulated image is the transformation of a sinus iamges I(x,y) =sin(x)
// by a mapping  X,Y  ->  (X + aY^2 , Y), we add also an attenaution functin
// that make image darker

void cAppliNewFrange::MakeImSimul()
{

    tREAL8 aDistIntrFr = 300.0 / mZoomRed; // Distance betweeb franges
    tREAL8 aMulY=10.0; // multipiler of the parabol X= MulY YN^2  with normalized YN
    tREAL8 aExp = 4.0;  // Exponent  of 1+sinus => the highest, give thinner franges
    tREAL8 aMiddleY = mSzRed.y() / 2.0; // Y of Middle horizonatl line

    for (const auto & aPix : *mDImRed)
    {
        tREAL8 aYNorm = (aPix.y()-aMiddleY) / aMiddleY;  // Nomalize Y : i.e. in [-1,1]

        // Phase  of sinus
        tREAL8 aPhase = aPix.x() -Square(aYNorm) * aMiddleY * aMulY;
        aPhase /= (aDistIntrFr /(2*M_PI));

        // compute peridic funtion in [0,1]
        tREAL8 aAmpl = std::max(0.0,(1+sin(aPhase)) /2.0);
        aAmpl = std::pow(aAmpl,aExp); // make thinner

        // make an attenuation , darker when we are further away of center
        aAmpl = aAmpl / (1 + std::pow(std::abs(aYNorm),2)*3.0 );

        mDImRed->SetV(aPix,aAmpl*255.0);  // now put it
    }
}


std::string cAppliNewFrange::NameVisu(const std::string & aPref) const
{
    return mPhProj.DirVisuAppli() + LastPrefix(mNameIm) + aPref + ".tif";
}


bool cAppliNewFrange::NewCC(bool isHoriz,const cPt2di& aPix)
{
    const std::vector<cPt2di> & a8Neigh =  Alloc8Neighbourhood();
    cDataIm2D<tU_INT1> & aImMax = isHoriz ? *mDImMaxHor : * mDImMaxVert;
    std::vector<cPt2di> aVCC;

    ConnectedComponent (aVCC,aImMax ,a8Neigh, aPix,1,2);

    tREAL8 aThrs = isHoriz ? 5 : 20 ;
    if (aVCC.size()<aThrs)
        return false;

    //cBox2di aBox(cTplBox<int,2>::FromVect(aVCC));
    cTplBox<int,2>  aBox =cTplBox<int,2>::FromVect(aVCC,true);

    bool aY0Inf =  (aBox.P0().y() <= mYC);
    bool aY1Sup =  (aBox.P1().y() > mYC);
    // MMVII_INTERNAL_ASSERT_always(aY0Inf<=aY1Sup," Ordeeerr in box");

     bool doCross = (aY0Inf && aY1Sup);
// FakeUseIt(doCross);
     cPt2di aDir(0,1);
     if (isHoriz)
     {
         // if dont cross
        if (!doCross)
            return false;
     }
     else
     {
      //   if (doCross)
       //     return false;
         aDir = (aBox.P0().y()<mYC)  ? cPt2di(1,0) : cPt2di(-1,0);
     }


     cWhichMinMax<cPt2di,tREAL8>  aWMM;
     cWeightAv<tREAL8,tREAL8> aWScore;

     for (const auto & aPt : aVCC)
     {
         aWMM.Add(aPt,Scal(aDir,aPt));

         tREAL8 aV = mDImRedBlur->GetV(aPix);
         tREAL8 aAmpl = mDImRadFrange->GetV(aPix.y()) -mRadiomBackGround;
         aV = (aV-mRadiomBackGround)/aAmpl;
         aV = std::clamp(aV,0.0,1.0);
         aWScore.Add(1.0,aV);
     }

     tREAL8 aScore = aWScore.Average();
     tREAL8 aScoreMixte = aScore * aVCC.size();
     bool isOk =    (aScore > 0.5)
                 || (aScoreMixte > 15)
                 || ((!isHoriz) && (aVCC.size() > 150))
            ;
     tSeg2dr aSeg(ToR(aWMM.IndMin()),ToR(aWMM.IndMax()));


     cConnComp aCC(aVCC,aSeg,isHoriz,aScore,isOk);
   //  cConnComp(const tSeg2dr & aSeg,bool isHor) :

     mListCC.push_back(aCC);

     return true;
}

void cAppliNewFrange::ConnecCompMaxLoc(bool isHoriz)
{
    cDataIm2D<tU_INT1> & aImMax = isHoriz ? *mDImMaxHor : * mDImMaxVert;

    for (const auto aPix : aImMax)
    {
        if (aImMax.GetV(aPix)==1)
        {
            NewCC(isHoriz,aPix);
        }
    }
}


bool cAppliNewFrange::IsMax(cPt2di aP0,cPt2di aDp0,int aNb) const
{
    tREAL4 aV0 = mDImRedBlur->GetV(aP0);
    for (int aK=1 ; aK<aNb ; aK++)
    {
        cPt2di aDP = aDp0 * aK;

        if (     ( mDImRedBlur->GetV(aP0 + aDP) > aV0)
              || ( mDImRedBlur->GetV(aP0 - aDP) > aV0)
           )
        {
            return false;
        }
    }

    return true;
}

void cAppliNewFrange::MakeImageMaxLoc()
{
    mImMaxHor = cIm2D<tU_INT1>(mSzRed,nullptr,eModeInitImage::eMIA_Null);
    mDImMaxHor = & (mImMaxHor.DIm());
    mImMaxVert = cIm2D<tU_INT1>(mSzRed,nullptr,eModeInitImage::eMIA_Null);
    mDImMaxVert = & (mImMaxVert.DIm());


    int aNbMaxX=8;
    int aNbMaxY=3;

    for (const auto aPix : mDImRed->Interior(1+std::max(aNbMaxX,aNbMaxY)))
    {
         if ( IsMax(aPix,cPt2di(1,0),aNbMaxX))
        {
             mDImMaxHor->SetV(aPix,1);
        }

        if (  IsMax(aPix,cPt2di(0,1),aNbMaxY))
        {
              mDImMaxVert->SetV(aPix,1);
        }
    }

    ConnecCompMaxLoc(true);
    ConnecCompMaxLoc(false);

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
    mImZ1 = tIm::FromFile(mNameIm);
    mDImZ1 = & (mImZ1.DIm());
    mImRed = mImZ1.BiCubicDeZoom(mZoomRed);
    mDImRed = &(mImRed.DIm());
    mSzRed = mDImRed->Sz();

    if (mDoSimul)
        MakeImSimul();


    mImRedBlur = mImRed.Dup();
    mDImRedBlur = & (mImRedBlur.DIm());
    ExpFilterOfStdDev(*mDImRedBlur,5,2.0,5.0);

    DoTensorProcessing();
    ComputeRadiomCste();

    MakeImageMaxLoc();



    if (mDoVisu)
    {
        tREAL8 aNbX=400;

        mImVisu = cRGBImage(mSzRed + cPt2di(aNbX,0),cRGBImage::White);
        cRGBImage aVisuMaxHor (mSzRed);
        cRGBImage aVisuMaxVert (mSzRed);
        cRGBImage aVisuArrow (mSzRed);

        cRGBImage aVisuTeta(mSzRed);
        for (const auto aPix : *mDImRed)
        {
            cPt2dr aTens(mDImTx->GetV(aPix),mDImTy->GetV(aPix));
            cPt2dr aRhoTeta = ToPolar(aTens,0.0);
            tREAL8 aTeta =  aRhoTeta.y();
            aVisuTeta.SetRGBPix(aPix,HSI_2_RGB(cPt3dr(aTeta,1.0,0.5)));
        }

        tElIm aVMin,mHighRadiom;
        GetBounds(aVMin,mHighRadiom,*mDImRed);

        for (const auto aPix : *mDImRed)
        {

            tREAL8 aRad = ((mDImRed->GetV(aPix)-mRadiomBackGround) /(mHighRadiom-mRadiomBackGround)) * 255.0;
            tINT4 aVal = std::clamp(round_ni(aRad),0,255);

            mImVisu.SetGrayPix(aPix,aVal);
            aVisuMaxHor.SetGrayPix(aPix,aVal);
            aVisuMaxVert.SetGrayPix(aPix,aVal);
            aVisuArrow.SetGrayPix(aPix,aVal);
        }






        for (const auto aPix : * mDImRed)
        {
            if ( mDImMaxLoc->GetV(aPix))
                aVisuMaxHor.SetRGBPix(aPix,cRGBImage::Yellow);
            if ( mDImMaxVert->GetV(aPix))
                aVisuMaxVert.SetRGBPix(aPix,cRGBImage::Cyan);
        }
        //aVisuMaxHor.DrawLine(cPt2dr(0,mYC),cPt2dr(mSzRed.x(),mYC),cRGBImage::Blue,1.0);
        aVisuMaxVert.DrawLine(cPt2dr(0,mYC),cPt2dr(mSzRed.x(),mYC),cRGBImage::Blue,1.0);
        aVisuArrow.DrawLine(cPt2dr(0,mYC),cPt2dr(mSzRed.x(),mYC),cRGBImage::Blue,1.0);



        mImVisu.DrawLine(cPt2dr(0,mYC),cPt2dr(mSzRed.x()+aNbX,mYC),cRGBImage::Blue,1.0);

        std::vector<cPt2dr> aVPtsIntegral;
        for(const auto aPtY : mImTgt.DIm())
        {
            int anY = aPtY.x();
            tREAL8 aXInt = mDImIntegr->GetV(anY);
            aVPtsIntegral.push_back(cPt2dr(aXInt+200.0,anY));
        }
        for(const auto aPtY : mImTgt.DIm())
        {
           int aXMil = mSzRed.x()+aNbX/2;
           int anY = aPtY.x();
           // show middel line
           mImVisu.SetRGBPix(cPt2di(aXMil,anY),cRGBImage::Green);


           cPt2di  aPtRad(round_ni(mSzRed.x()+mDImRadFrange->GetV(anY)*0.5),anY);
           mImVisu.SetRGBPix(aPtRad,cRGBImage::Gray128);


           cPt2di  aPtTgt(round_ni(aXMil+mDImTgt->GetV(anY)*10.0),anY);
           mImVisu.SetRGBPix(aPtTgt,cRGBImage::Red);
           cPt2di  aPtTeta(round_ni(aXMil+mDImTeta->GetV(anY)*25.0),anY);
           mImVisu.SetRGBPix(aPtTeta,cRGBImage::Blue);



           // Show image integrale in image
           //int aXInt = mDImIntegr->GetV(anY);
           //mImVisu.SetRGBPix(cPt2di(aXInt+100,anY),cRGBImage::Red);
           if (anY>0)
           {
               cPt2dr aP1 = aVPtsIntegral.at(anY) ;
               cPt2dr aP2 = aVPtsIntegral.at(anY-1);
               if (mImVisu.InsideBL(aP1) && mImVisu.InsideBL(aP2))
               {
                 // StdOut() << "PTTTT " << aP1 << aP2 << "\n";
                  mImVisu.DrawLine(aP1,aP2,cRGBImage::Red);
               }
           }
        }


        // "ARROW"  visu
        for (const auto & aCC : mListCC)
        {
            const tSeg2dr& aSeg = aCC.mSeg;
            cPt3di aCol = aCC.mIsHor ? cRGBImage::Yellow : cRGBImage::Cyan;
            if (!aCC.mIsOk)
                aCol = cRGBImage::Magenta;
            if (true) // (aCC.mIsOk)
            {
                for (const auto & aPix : aCC.mPts)
                    aVisuArrow.SetRGBPix(aPix,aCol);
                aVisuArrow.DrawCircle(cRGBImage::Red,aSeg.P1(),3.0);
                aVisuArrow.DrawCircle(cRGBImage::Green,aSeg.P2(),3.0);
            }
        }



        if (mDoVisu>=2)
        {
           aVisuMaxVert.ToFile(NameVisu("ImMaxVert"));
           mDImRedBlur->ToFile(NameVisu("Blured"));
           aVisuTeta.ToFile(NameVisu("TetaTens"));
           aVisuMaxHor.ToFile(NameVisu("ImMaxHor"));
           aVisuArrow.ToFile(NameVisu("ImArrow"));
        }
        mImVisu.ToFile(NameVisu("ImRed"));


     }
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
