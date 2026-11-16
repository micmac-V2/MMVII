//#include "MMVII_PCSens.h"

// Test commit

#include "cAppliNewFranges.h"
#include "MMVII_Tpl_Images.h"
#include "MMVII_Tpl_NormImage.h"

namespace MMVII
{

namespace NS_FrangesDetect
{


/* =================================================== */
/*              Specific functions                     */
/* =================================================== */



std::string cAppliNewFrange::NameVisu(const std::string & aPref) const
{
    return mPhProj.DirVisuAppli() + LastPrefix(mNameIm) + aPref + ".tif";
}




// The simulated image is the transformation of a sinus iamges I(x,y) =sin(x)
// by a mapping  X,Y  ->  (X + aY^2 , Y), we add also an attenaution functin
// that make image darker

void cAppliNewFrange::MakeImSimul()
{

    tREAL8 aDistIntrFr = 300.0 / mZoomRed; // Distance betweeb franges
    tREAL8 aMulY=10.0; // multipiler of the parabol X= MulY YN^2  with normalized YN
    tREAL8 aExp = 4.0;  // Exponent  of 1+sinus => the highest, give thinner franges

    if (mDoSimul==1)
    {
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
    else if (mDoSimul==2)
    {
        for (const auto & aPix : *mDImRed)
        {
            cPt2dr aRhoTeta = ToPolar(ToR(aPix) - ToR(mSzRed)/2.0,0.0);  // Nomalize Y : i.e. in [-1,1]

            // Phase  of sinus
            tREAL8 aPhase = aRhoTeta.x();
            aPhase /= (aDistIntrFr /(2*M_PI));

            // compute peridic funtion in [0,1]
            tREAL8 aAmpl = std::max(0.0,(1+sin(aPhase)) /2.0);
            aAmpl = std::pow(aAmpl,aExp); // make thinner

            mDImRed->SetV(aPix,aAmpl*255.0);  // now put it
        }
    }
}

void cAppliNewFrange::ShowCurvesIntegral(cRGBImage& anIm,bool doShowC,bool doShowSep)
{
    for (const auto aML : mVMaxLoc)
    {
        for(const auto aPtY : mImTgt.DIm())
        {
            int anY = aPtY.x();
            if (anY>0)
            {
                tREAL8 aDX0 = mDImIntegr->GetV(anY-1);
                tREAL8 aDX1 = mDImIntegr->GetV(anY);


                if  (doShowC && (aDX1+aML.mXC<mSzRed.x()))
                    anIm.DrawLine(cPt2dr(aDX0+aML.mXC,anY-1),cPt2dr(aDX1+aML.mXC,anY),cRGBImage::Red);

                if (doShowSep && (aDX1+aML.mX0<mSzRed.x()))
                    anIm.DrawLine(cPt2dr(aDX0+aML.mX0,anY-1),cPt2dr(aDX1+aML.mX0,anY),cRGBImage::Blue);
                if (doShowSep && (aDX1+aML.mX1<mSzRed.x()))
                    anIm.DrawLine(cPt2dr(aDX0+aML.mX1,anY-1),cPt2dr(aDX1+aML.mX1,anY),cRGBImage::Blue);
             }
        }
    }
}


const cPt3di *cAppliNewFrange::CoulOfLabel(const cPt2di & aPix) const
{
    if ( LabelIsOn(aPix,eTypePix::eMaxTens))
        return  LabelIsOn(aPix,eTypePix::eProlong) ? &cRGBImage::Red : &cRGBImage::Green;

    if ( LabelIsOn(aPix,eTypePix::eMinTens))
        return  LabelIsOn(aPix,eTypePix::eProlong) ? &cRGBImage::Orange : &cRGBImage::Blue;

    return nullptr;

    if ( LabelIsOn(aPix,eTypePix::eTypeSmalCC))
        return & cRGBImage::Yellow;

    if ( LabelIsOn(aPix,eTypePix::eTypeBadAvgCC))
        return & cRGBImage::Cyan;

    return nullptr;
}

cRGBImage cAppliNewFrange::ImBckGround(tREAL8 aMul) const
{
    return cRGBImage( Im2Make8BByGaussModel(aMul,*mDImRedBlur).DIm());
}

void cAppliNewFrange::SetImageLabel(cRGBImage & anIm,tREAL8 aCoefAB) const
{
    for (const auto aPix :  mDImRegion->Interior(1))
    {
           int aN1 = mDImRegion->GetV(aPix);
           if (aN1>=0)
           {
               int aNum1CC0 = mGrReg.VertexOfNum(aN1).Attr().mNumCC0;

               anIm.SetRGBPix(aPix,cRGBImage::KthRandomCol(mImGrayBG.DIm().GetV(aPix),aNum1CC0,0.3));
           }

    }
}

void cAppliNewFrange::ShowExtres(cRGBImage & anIm) const
{
    for (const auto & aV : mGrReg.AllVertices())
    {
        const auto & anAttr = aV->Attr();
       /* for (const auto & aPix : anAttr.mVPts)
            anIm.SetRGBPix(aPix,cRGBImage::Red);*/

        int aNL = aV->Attr().mNumMaxLoc;
        bool hasNum = aNL >= 0;

        // show pathes
        cPt3di aColPath =  cRGBImage::Green;
        if (hasNum)
        {
            if (aV->Attr().mProvMaxLoc==0)
            {
                aColPath = cRGBImage::Red;
            }
            else if (aV->Attr().mProvMaxLoc==1)
            {
                aColPath = cRGBImage::Orange;
            }
            else if (aV->Attr().mProvMaxLoc==2)
            {
                aColPath = cRGBImage::Yellow;
            }
        }
        for (const auto & aPix : anAttr.mPath)
            anIm.SetRGBPix(aPix,aColPath);

        // show approx poly
        for (const auto anInd: anAttr.mIndApP)
        {
            anIm.DrawCircle(cRGBImage::Cyan,ToR(anAttr.mPath.at(anInd)),2.0);
        }


        // show extrmeties with 2 diff colors
        for (int aK=0 ; aK<2 ; aK++)
        {
            cPt2di aPt  = (aK==0) ? anAttr.mP1 : anAttr.mP2 ;
            cPt3di aCol = (aK==0) ? cRGBImage::Green :  cRGBImage::Red;

            for (const auto aRay : {3.0} ) //  {2.0,4.0,6.0})
            {
                  anIm.DrawCircle(aCol,ToR(aPt),aRay);
            }
        }

        if (hasNum)
        {
            anIm.DrawString(ToStr(aNL),cRGBImage::Orange,ToR(aV->Attr().mCenter),cPt2dr(0,0));
        }


        //  anIm.DrawCircle(cRGBImage::Orange,ToR(anAttr.mCenter),4.0);

    }

}


void cAppliNewFrange::DoVisu()
{
    if (!IsInit(&mPatVisu))
        return;

    mImGrayBG = Im2Make8BByGaussModel(0.5,*mDImRedBlur);
    mImRgbBG = cRGBImage(mImGrayBG.DIm());



    //  ------------ Generate image of Label by colour of reg -----------------------
    static const std::string TheNameLabel = "Label";
    // if (MatchRegex(TheNameLabel,mPatVisu))
    {
       cRGBImage anImLabel =    mImRgbBG.Dup();
       SetImageLabel(anImLabel,0.3);
       GenVisu(anImLabel,TheNameLabel);
    }

    //   ------------- Image of extre +
    {
         cRGBImage aVisuOrientedLines = mImRgbBG.Dup();

         ShowExtres(aVisuOrientedLines);
         ShowCurvesIntegral(aVisuOrientedLines,false,true);
         GenVisu(aVisuOrientedLines,"Extres");
    }





    // ------------------- Generate visualisation of tensor -----------------------
    static const std::string TheNameTetaTens = "TetaTens";
    static const std::string TheNameNormTens = "NormTens";

    //if (MatchRegex(TheNameLabel,mPatVisu))
    {
        cIm2D<tREAL8> aINorm(mSzRed);
        cRGBImage aVisuTeta(mSzRed);
        for (const auto aPix : *mDImRed)
        {
            cPt2dr aTens = mImTens.GradR(aPix);
            cPt2dr aRhoTeta = ToPolar(aTens,0.0);
            tREAL8 aTeta =  aRhoTeta.y();
            aVisuTeta.SetRGBPix(aPix,HSI_2_RGB(cPt3dr(aTeta,1.0,0.5)));
            aINorm.DIm().SetV(aPix,aRhoTeta.x());
        }
        GenVisu(aVisuTeta,TheNameTetaTens);
        GenVisu(Im2Make8BByMinMax(aINorm.DIm()).DIm(),TheNameNormTens);
    }





    //  Generate the back gound with central line & lims




    {
        cRGBImage aVisuCC0 =   mImRgbBG.Dup();

       // aVisuMTReg.DrawLine(cPt2dr(0,mYC),cPt2dr(mSzRed.x(),mYC),cRGBImage::Gray128,1.0);
        aVisuCC0.DrawLine(cPt2dr(0,mYC),cPt2dr(mSzRed.x(),mYC),cRGBImage::Gray128,1.0);


        for (const auto aPix :  *mDImRegion)
        {
            // show crest-line in green
           if (LabelIsOn(aPix,eTypePix::eMaxTens))
           {
            //   aVisuMTReg.SetRGBPix(aPix,cRGBImage::Green);
               aVisuCC0.SetRGBPix(aPix,cRGBImage::Green);
           }
           // show the frontier in cyan
           for (int aK=0 ; aK<2 ; aK++)
           {
               int aN1 = mDImRegion->GetV(aPix);
               cPt2di aPix2 = (aPix+FreemanV4[aK]);
               if (mDImRegion->Inside(aPix2))
               {
                   int aN2 = mDImRegion->GetV(aPix+FreemanV4[aK]);

                   if ((aN1!=aN2) && (aN1>=0) && (aN2>=0))
                   {
                     //  aVisuMTReg.SetRGBPix(aPix,cRGBImage::Cyan);
                       tVertexReg & aVert1 = mGrReg.VertexOfNum(aN1);
                       tVertexReg & aVert2 = mGrReg.VertexOfNum(aN2);

                       int aNum1CC0 = aVert1.Attr().mNumCC0;
                       int aNum2CC0 = aVert2.Attr().mNumCC0;

                       cPt3di aCol = cRGBImage::Red;
                       if (aNum1CC0!=aNum2CC0)
                       {
                           aVisuCC0.SetRGBPix(aPix,aCol);

                       }
                       else
                       {
                           /*
                           tEdgeReg * anE =aVert1.EdgeOfSucc(aVert2,true);
                           aCol = (anE && anE->AttrSym().mInMST) ? cRGBImage::Blue : cRGBImage::Cyan ;
                           */
                       }
                   }
               }
           }
           // show the valley in blue (partially overwrite valley)
           //if (LabelIsOn(aPix,eTypePix::eMinTens))
           //    aVisuMTReg.SetRGBPix(aPix,cRGBImage::Blue);
        }

        // show connexion
        /*for (const auto & anEdge : mGrReg.AllEdges_DirInit())
        {
             cPt2dr aP1 = ToR(anEdge->AttrOriented().mPt);
             cPt2dr aP2 = ToR(anEdge->AttrInvOr().mPt);
             cPt2dr aMil = (aP1+aP2) / 2.0;
             cPt2dr aV12 = aP2-aP1;
             tREAL8 aL = 5.0;
             cPt3di aColor = anEdge->AttrSym().mInMST ? cRGBImage::Red : cRGBImage::Yellow;
             aVisuMTReg.DrawLine(aMil-aV12*aL,aMil+aV12*aL,aColor);
        }*/

      /*  if (0)
        {
            for (const auto & aV : mGrReg.AllVertices())
            {
                const auto & a2E = aV->Attr().mExtres;
                for (const auto aRay : {3.0,5.0,7.0})
                {
                    aVisuCC0.DrawCircle(cRGBImage::Magenta,ToR(a2E.first) ,aRay);
                    aVisuCC0.DrawCircle(cRGBImage::Magenta,ToR(a2E.second),aRay);
                }
            }
        }


        for (const auto aML : mVMaxLoc)
        {
            aVisuMTReg.DrawCircle(cRGBImage::Red,cPt2dr(aML.mXC,mYC),3.0);
        }*/
        //GenVisu(aVisuMTReg,"MaxTReg");
        GenVisu(aVisuCC0,"CC0");

    }


    mImRgbBG.DrawLine(cPt2dr(0,mMinZYC),cPt2dr(mSzRed.x(),mMinZYC),cRGBImage::Green,1.0);
    mImRgbBG.DrawLine(cPt2dr(0,mMaxZYC),cPt2dr(mSzRed.x(),mMaxZYC),cRGBImage::Green,1.0);
    mImRgbBG.DrawLine(cPt2dr(0,mYC),cPt2dr(mSzRed.x(),mYC),cRGBImage::Red,1.0);

    GenVisu(mImRgbBG,"BackGround");


#if (0)
    if (0)
    {
        tIm aImFiltr = mImRed.Dup();
        ExpFilterOfStdDev(aImFiltr.DIm(),1,mIntervMax/10.0,1.5);

        cIm1D<tREAL4> aImL(mSzRed.x());
        for (int anY =0 ; anY<mSzRed.y() ; anY++)
        {
            aImFiltr.DIm().WriteLine(aImL.DIm(),anY);
           // NormalizeByHistoEq(aImL.DIm(),1.0,aImL.DIm(),mSzRed.x());
            FixDynByGaussModel(aImL.DIm(),0,1,aImL.DIm());
            aImFiltr.DIm().ReadLine(aImL.DIm(),anY);

        }
        //template<class T> cIm2D<tU_INT1> Im2Make8BByGaussModel(const cDataIm2D<T> & aIm,tREAL8 aEpsilon=1e-5)

         cIm2D<tU_INT1> aIm8B = Im2Make8BByMinMax(aImFiltr.DIm());
        GenVisu(aIm8B.DIm(),"FiltX");

        //MakeGrapheTensor(aIm8B);
    }

    // Generate images of angle of tensor using Hue (teinte)
    {
        cIm2D<tREAL8> aINorm(mSzRed);
        cRGBImage aVisuTeta(mSzRed);
        for (const auto aPix : *mDImRed)
        {
            //cPt2dr aTens(mDImTx->GetV(aPix),mDImTy->GetV(aPix));
            cPt2dr aTens = mImTens.GradR(aPix);
            cPt2dr aRhoTeta = ToPolar(aTens,0.0);
            tREAL8 aTeta =  aRhoTeta.y();
            aVisuTeta.SetRGBPix(aPix,HSI_2_RGB(cPt3dr(aTeta,1.0,0.5)));
            aINorm.DIm().SetV(aPix,aRhoTeta.x());
        }
        GenVisu(aVisuTeta,"TetaTens");
        GenVisu(Im2Make8BByMinMax(aINorm.DIm()).DIm(),"NormTens");

    }

    tREAL8 aRabXDetail=400;
    tREAL8 aRabYDetail=200;


    //  Image with several detail : tang/ teta-tens / integrale / Rad(Y) of Frange
    cPt2di aSzImDetail = mSzRed + cPt2di(aRabXDetail,aRabYDetail);
    cRGBImage anImDetail = cRGBImage(aSzImDetail,cRGBImage::White);
    anImDetail.FillRectangle(cRGBImage::Yellow,mSzRed,aSzImDetail,cPt3dr::PCste(0.0));

    // Image with normalzized dyn that will be used as background
    cRGBImage aImGrayBackGround (mSzRed);

    // image that visualize oriented chain (+- arrow ...)
    cRGBImage aVisuArrow (mSzRed);
    cRGBImage aVisuArrowUnif (mSzRed);


    cIm2D<tU_INT1> aImNorm = Im2Make8BByGaussModel(*mDImRedBlur);
    // initialize the raidometry with some enhancement
    for (const auto aPix : *mDImRed)
    {
        tREAL8 aHR = mDImRadFrange->GetV(aPix.y());
        tREAL8 aRad = ((mDImRedBlur->GetV(aPix)-mRadiomBackGround) /(aHR-mRadiomBackGround)) * 255.0;
        tINT4 aVal = std::clamp(round_ni(aRad),0,255);

        aVal = aImNorm.DIm().GetV(aPix);
        aImGrayBackGround.SetGrayPix(aPix,aVal);

        anImDetail.SetGrayPix(aPix,aVal);
        aVisuArrow.SetGrayPix(aPix,aVal);
        aVisuArrowUnif.SetRGBPix(aPix,cRGBImage::Gray128);
    }


    {
        cRGBImage aVisuMaxTens = aImGrayBackGround.Dup();
        aVisuMaxTens.DrawLine(cPt2dr(0,mYC),cPt2dr(mSzRed.x(),mYC),cRGBImage::Gray128,1.0);
        for (const auto aPix : * mDImRed)
        {
            const cPt3di * aCoul = CoulOfLabel(aPix) ;
            if (aCoul)
                aVisuMaxTens.SetRGBPix(aPix,*aCoul);
        }
        GenVisu(aVisuMaxTens,"MaxTens");
    }

    {
        cRGBImage aImReg (mSzRed);
        for (const auto aPix : * mDImRegion)
        {
            int aNum = mDImRegion->GetV(aPix);
            if (aNum>=0)
            {
               // const auto & anAttr = mGrReg.VertexOfNum(aNum).Attr();
                aImReg.SetRGBPix(aPix,cRGBImage::KthRandomCol(aNum));
            }
        }
        GenVisu(aImReg,"Region");
    }
    {
        cRGBImage aVisuCC0 =   aImGrayBackGround.Dup();


        aVisuMTReg.DrawLine(cPt2dr(0,mYC),cPt2dr(mSzRed.x(),mYC),cRGBImage::Gray128,1.0);
        aVisuCC0.DrawLine(cPt2dr(0,mYC),cPt2dr(mSzRed.x(),mYC),cRGBImage::Gray128,1.0);


        for (const auto aPix :  *mDImRegion)
        {
            // show crest-line in green
           if (LabelIsOn(aPix,eTypePix::eMaxTens))
           {
               aVisuMTReg.SetRGBPix(aPix,cRGBImage::Green);
               aVisuCC0.SetRGBPix(aPix,cRGBImage::Green);
           }
           // show the frontier in cyan
           for (int aK=0 ; aK<2 ; aK++)
           {
               int aN1 = mDImRegion->GetV(aPix);
               cPt2di aPix2 = (aPix+FreemanV4[aK]);
               if (mDImRegion->Inside(aPix2))
               {
                   int aN2 = mDImRegion->GetV(aPix+FreemanV4[aK]);


                   if ((aN1!=aN2) && (aN1>=0) && (aN2>=0))
                   {
                       aVisuMTReg.SetRGBPix(aPix,cRGBImage::Cyan);
                       tVertexReg & aVert1 = mGrReg.VertexOfNum(aN1);
                       tVertexReg & aVert2 = mGrReg.VertexOfNum(aN2);

                       int aNum1CC0 = aVert1.Attr().mNumCC0;
                       int aNum2CC0 = aVert2.Attr().mNumCC0;

                       cPt3di aCol = cRGBImage::Red;
                       if (aNum1CC0!=aNum2CC0)
                       {

                       }
                       else
                       {
                           tEdgeReg * anE =aVert1.EdgeOfSucc(aVert2,true);
                           aCol = (anE && anE->AttrSym().mInMST) ? cRGBImage::Blue : cRGBImage::Cyan ;

                       }
                       aVisuCC0.SetRGBPix(aPix,aCol);
                   }
               }
           }
           // show the valley in blue (partially overwrite valley)
           if (LabelIsOn(aPix,eTypePix::eMinTens))
               aVisuMTReg.SetRGBPix(aPix,cRGBImage::Blue);
        }

        // show connexion
        for (const auto & anEdge : mGrReg.AllEdges_DirInit())
        {
             cPt2dr aP1 = ToR(anEdge->AttrOriented().mPt);
             cPt2dr aP2 = ToR(anEdge->AttrInvOr().mPt);
             cPt2dr aMil = (aP1+aP2) / 2.0;
             cPt2dr aV12 = aP2-aP1;
             tREAL8 aL = 5.0;
             cPt3di aColor = anEdge->AttrSym().mInMST ? cRGBImage::Red : cRGBImage::Yellow;
             aVisuMTReg.DrawLine(aMil-aV12*aL,aMil+aV12*aL,aColor);
        }

        if (0)
        {
            for (const auto & aV : mGrReg.AllVertices())
            {
                const auto & a2E = aV->Attr().mExtres;
                for (const auto aRay : {3.0,5.0,7.0})
                {
                    aVisuCC0.DrawCircle(cRGBImage::Magenta,ToR(a2E.first) ,aRay);
                    aVisuCC0.DrawCircle(cRGBImage::Magenta,ToR(a2E.second),aRay);
                }
            }
        }


        for (const auto aML : mVMaxLoc)
        {
            aVisuMTReg.DrawCircle(cRGBImage::Red,cPt2dr(aML.mXC,mYC),3.0);
        }
        GenVisu(aVisuMTReg,"MaxTReg");
        GenVisu(aVisuCC0,"CC0");

    }

    GenVisu(*mImGrad.mDGx,"GradX");
    GenVisu(*mImGrad.mDGy,"GradY");

    //  Visualize the central line
    anImDetail.DrawLine(cPt2dr(0,mYC),cPt2dr(mSzRed.x()+aRabXDetail,mYC),cRGBImage::Cyan,1.0);



    // Transformate vector Integal in points, arbirary origin


//    template <class Type> std::pair<cIm1D<Type>,cIm1D<Type>> ComputeAvgAndDev(const cDataIm1D<Type> & aDI1,int  aNbIt,double aStdDev);

    tIm1D aImVisuNT = Im1DFixDynByMax(*mDImCumNormT,100.0);
    tREAL8 aSigma = 2.0;
    auto [aImMoy,aImDev] = ComputeAvgAndDev(*mDImTeta,3,aSigma);

    //  ---------- Make image details ----------------------
    for(const auto aPtY : mImTgt.DIm())
    {
       int aXMil = mSzRed.x()+aRabXDetail/2;
       int anY = aPtY.x();

       // ----------------- Curves in the white rectangle ------------------------

           //   ---  show middle line  ---
       anImDetail.SetRGBPix(cPt2di(aXMil,anY),cRGBImage::Green);

           //   ---  show curves that estimate radiomtry of franges ------------
       cPt2di  aPtRad(round_ni(mSzRed.x()+mDImRadFrange->GetV(anY)*0.5),anY);
       anImDetail.SetRGBPix(aPtRad,cRGBImage::Gray128);

           // --- show image of norm of tensor ---
       if (0)
       {
          cPt2di  aPNormFr(round_ni(mSzRed.x()+aImVisuNT.DIm().GetV(anY)+2),anY);
          anImDetail.SetRGBPix(aPNormFr,cRGBImage::Green);
       }
       {
           cPt2di  aPNormFr(round_ni(mSzRed.x()+aImDev.DIm().GetV(anY)*500.0/aSigma+2),anY);
           anImDetail.SetRGBPix(aPNormFr,cRGBImage::Green);
       }

           //   ---  curves of teta and tangent --------
       cPt2di  aPtTgt(round_ni(aXMil+mDImTgt->GetV(anY)*10.0),anY);
       anImDetail.SetRGBPix(aPtTgt,cRGBImage::Red);
       cPt2di  aPtTeta(round_ni(aXMil+mDImTeta->GetV(anY)*25.0),anY);
       anImDetail.SetRGBPix(aPtTeta,cRGBImage::Blue);

    }
    for(const auto aPtX : *mDImProjInt)
    {
        int anX = aPtX.x();
        tREAL8 anY  = ((mDImProjInt->GetV(anX)-mMinDPI) / (mMaxDPI-mMinDPI)) * aRabYDetail;
        anImDetail.SetRGBPix(cPt2di(anX,mSzRed.y() + anY *0.95 + 5),cRGBImage::Black);
    }

   ShowCurvesIntegral(anImDetail);



    //  ------------  visualisation of oriented chains -----------------------

   /*
    for (auto aIm : {aVisuArrow,aVisuArrowUnif})
    {
        aIm.DrawLine(cPt2dr(0,mYC),cPt2dr(mSzRed.x(),mYC),cRGBImage::Blue,1.0);
        for (const auto & aCC : mVectCC)
        {
            const tSeg2dr& aSeg = aCC->mSeg;
            cPt3di aCol = aCC->mIsHor ? cRGBImage::Yellow : cRGBImage::Cyan;
            if (!aCC->mIsOk)
                aCol = cRGBImage::Magenta;
            if (aCC->mIsOk)
            {
                for (const auto & aPix : aCC->mPts)
                    aIm.SetRGBPix(aPix,aCol);
                aIm.DrawCircle(cRGBImage::Red,aSeg.P1(),3.0);
                aIm.DrawLine(aSeg.P1(),aSeg.P1()+aCC->mTg1*10.0,cRGBImage::Red);
                aIm.DrawCircle(cRGBImage::Green,aSeg.P2(),3.0);
                aIm.DrawLine(aSeg.P2(),aSeg.P2()+aCC->mTg2*10.0,cRGBImage::Green);

            }
        }
    }*/


    GenVisu(anImDetail,"Details");
    GenVisu(aVisuArrow,"ArrowIm");
    GenVisu(aVisuArrowUnif,"ArrowUnif");

    GenVisu(*mDImRedBlur,"Blured");


    GenVisu(aImNorm.DIm(),"Normed");

#endif

    // Possibly generates names of visu (if PatInit & NoVisu Gene & 1 single file)

    if (IsInit(&mPatVisu) && ((mNbVisuGen==0)||(mPatVisu==".*")) && (LevelCall()==0)  )
    {
        StdOut()  <<  "--OPTION VISU=\n" ;
        for (const auto & aNameV : mVNameVisu)
            StdOut() << "  * " << aNameV<< "\n";
    }
}

};  // NS_FrangesDetect



}; //  namespace MMVII
