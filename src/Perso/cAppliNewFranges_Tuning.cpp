//#include "MMVII_PCSens.h"

// Test commit

#include "cAppliNewFranges.h"




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


void cAppliNewFrange::DoVisu()
{
    if (!IsInit(&mPatVisu))
        return;


    // Generate images of angle of tensor using Hue (teinte)
    {
        cRGBImage aVisuTeta(mSzRed);
        for (const auto aPix : *mDImRed)
        {
            //cPt2dr aTens(mDImTx->GetV(aPix),mDImTy->GetV(aPix));
            cPt2dr aTens = mImTens.GradR(aPix);
            cPt2dr aRhoTeta = ToPolar(aTens,0.0);
            tREAL8 aTeta =  aRhoTeta.y();
            aVisuTeta.SetRGBPix(aPix,HSI_2_RGB(cPt3dr(aTeta,1.0,0.5)));
        }
        GenVisu(aVisuTeta,"TetaTens");
    }

    tREAL8 aNbX=400;

    //  Image with several detail : tang/ teta-tens / integrale / Rad(Y) of Frange
    cRGBImage anImDetail = cRGBImage(mSzRed + cPt2di(aNbX,0),cRGBImage::White);
    // image that visualize oriented chain (+- arrow ...)
    cRGBImage aVisuArrow (mSzRed);
    cRGBImage aVisuArrowUnif (mSzRed);



    // initialize the raidometry with some enhancement

    for (const auto aPix : *mDImRed)
    {
        tREAL8 aHR = mDImRadFrange->GetV(aPix.y());
        tREAL8 aRad = ((mDImRed->GetV(aPix)-mRadiomBackGround) /(aHR-mRadiomBackGround)) * 255.0;
        tINT4 aVal = std::clamp(round_ni(aRad),0,255);

        anImDetail.SetGrayPix(aPix,aVal);
        aVisuArrow.SetGrayPix(aPix,aVal);
        aVisuArrowUnif.SetRGBPix(aPix,cRGBImage::Gray128);

    }


    //  generate images of max loc  Horiz/Vert/Tens and Mixte
    {
        cRGBImage aVisuMaxHor = aVisuArrow.Dup();
        cRGBImage aVisuMaxVert = aVisuArrow.Dup();
        cRGBImage aVisuMaxTens = aVisuArrow.Dup();
        cRGBImage aVisuMaxMixte = aVisuArrow.Dup();

        for (const auto aPix : * mDImRed)
        {
            if ( mDImMaxLocTD->GetV(aPix))
            {
                aVisuMaxTens.SetRGBPix(aPix,cRGBImage::Green);
                aVisuMaxMixte.SetRGBPix(aPix,cRGBImage::Green);
            }
            if ( mDImMaxHor->GetV(aPix))
            {
                aVisuMaxHor.SetRGBPix(aPix,cRGBImage::Yellow);
                aVisuMaxMixte.SetRGBPix(aPix,cRGBImage::Yellow);
            }
            if ( mDImMaxVert->GetV(aPix))
            {
                aVisuMaxVert.SetRGBPix(aPix,cRGBImage::Cyan);
                aVisuMaxMixte.SetRGBPix(aPix,cRGBImage::Cyan);
            }
        }


        GenVisu(aVisuMaxTens,"MaxLocTens");
        GenVisu(aVisuMaxHor,"MaxLocHor");
        GenVisu(aVisuMaxVert,"MaxLocVert");
        GenVisu(aVisuMaxMixte,"MaxLocMixte");
    }

    //  Visualize the central line
    anImDetail.DrawLine(cPt2dr(0,mYC),cPt2dr(mSzRed.x()+aNbX,mYC),cRGBImage::Blue,1.0);



    // Transformate vector Integal in points, arbirary origin
    std::vector<cPt2dr> aVPtsIntegral;
    for(const auto aPtY : mImTgt.DIm())
    {
        int anY = aPtY.x();
        tREAL8 aXInt = mDImIntegr->GetV(anY);
        aVPtsIntegral.push_back(cPt2dr(aXInt+200.0,anY));
    }


    //  ---------- Make image details ----------------------
    for(const auto aPtY : mImTgt.DIm())
    {
       int aXMil = mSzRed.x()+aNbX/2;
       int anY = aPtY.x();

       // ----------------- Curves in the white rectangle ------------------------

           //   ---  show middle line  ---
       anImDetail.SetRGBPix(cPt2di(aXMil,anY),cRGBImage::Green);

           //   ---  show curves that estimate radiomtry of franges ------------
       cPt2di  aPtRad(round_ni(mSzRed.x()+mDImRadFrange->GetV(anY)*0.5),anY);
       anImDetail.SetRGBPix(aPtRad,cRGBImage::Gray128);

           //   ---  curves of teta and tangent --------
       cPt2di  aPtTgt(round_ni(aXMil+mDImTgt->GetV(anY)*10.0),anY);
       anImDetail.SetRGBPix(aPtTgt,cRGBImage::Red);
       cPt2di  aPtTeta(round_ni(aXMil+mDImTeta->GetV(anY)*25.0),anY);
       anImDetail.SetRGBPix(aPtTeta,cRGBImage::Blue);

       // ------------------- Show image integrale in image ---------------------
       if (anY>0)
       {
           cPt2dr aP1 = aVPtsIntegral.at(anY) ;
           cPt2dr aP2 = aVPtsIntegral.at(anY-1);
           if (mDImRedBlur->InsideBL(aP1) && mDImRedBlur->InsideBL(aP2))
           {
              anImDetail.DrawLine(aP1,aP2,cRGBImage::Red);
           }
       }
    }


    //  ------------  visualisation of oriented chains -----------------------

    for (auto aIm : {aVisuArrow,aVisuArrowUnif})
    {
        aIm.DrawLine(cPt2dr(0,mYC),cPt2dr(mSzRed.x(),mYC),cRGBImage::Blue,1.0);
        for (const auto & aCC : mListCC)
        {
            const tSeg2dr& aSeg = aCC.mSeg;
            cPt3di aCol = aCC.mIsHor ? cRGBImage::Yellow : cRGBImage::Cyan;
            if (!aCC.mIsOk)
                aCol = cRGBImage::Magenta;
            if (aCC.mIsOk)
            {
                for (const auto & aPix : aCC.mPts)
                    aIm.SetRGBPix(aPix,aCol);
                aIm.DrawCircle(cRGBImage::Red,aSeg.P1(),3.0);
                aIm.DrawCircle(cRGBImage::Green,aSeg.P2(),3.0);
            }
        }
    }


    GenVisu(anImDetail,"Details");
    GenVisu(aVisuArrow,"ArrowIm");
    GenVisu(aVisuArrowUnif,"ArrowUnif");

    GenVisu(*mDImRedBlur,"Blured");

    // Possibly generates names of visu (if PatInit & NoVisu Gene & 1 single file)

    if (IsInit(&mPatVisu) && (mNbVisuGen==0) && (LevelCall()==0)  )
    {
        StdOut()  <<  "--OPTION VISU=\n" ;
        for (const auto & aNameV : mVNameVisu)
            StdOut() << "  * " << aNameV<< "\n";
    }
}

};  // NS_FrangesDetect



}; //  namespace MMVII
