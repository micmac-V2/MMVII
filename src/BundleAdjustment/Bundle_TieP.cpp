#include "BundleAdjustment.h"
#include "MMVII_Stringifier.h"


namespace MMVII
{

    /* ---------------------------------------- */
    /*            Tie Points                    */
    /* ---------------------------------------- */

cBA_TieP::cBA_TieP(const std::string & aName,cComputeMergeMulTieP* aMTP,const cStdWeighterResidual &aResW):
    mName           (aName),
    mMTP            (aMTP),
    mTieP_Weighter  (aResW)
{
}

cBA_TieP::~cBA_TieP()
{
    delete mMTP;
}


    /* ---------------------------------------- */
    /*            Tie Points                    */
    /* ---------------------------------------- */

void cMMVII_BundleAdj::AddMTieP(const std::string & aName,cComputeMergeMulTieP  * aMTP,const cStdWeighterResidual & aWIm)
{
    mVTieP.push_back(new cBA_TieP(aName,aMTP,aWIm));
}

void cMMVII_BundleAdj::SetTiePShowPerMil(const std::vector<int> & aTiePShowPerMil)
{
     mTiePShowPerMil = aTiePShowPerMil;
}



void cMMVII_BundleAdj::OneItere_TieP(const cBA_TieP& aBA_TieP)
{
   cComputeMergeMulTieP * aMTP = aBA_TieP.mMTP;
   const cStdWeighterResidual&    aTieP_Weighter = aBA_TieP.mTieP_Weighter;

   // update the bundle point by 3D-intersection:
   // To see : maybe don't update each time; probably add some robust option
   aMTP->SetPGround();

   cWeightAv<tREAL8> aWeigthedRes;
   cStdStatRes aStatResidual;
   bool doStatRes = ! mTiePShowPerMil.empty();

   cWeightAv<tREAL8> aUnWeightedRed; ///< MPD : a conserver, ils ont mysterieusement disparus; mais il sont pourtant tres utiles
   for (const auto & aPair : aMTP->Pts())
   {
       const auto & aConfig  = aPair.first;

       // local vector of sensor & colinearity equation, directly indexale in [0,NbIm]
       std::vector<cSensorImage *> aVS ;
       std::vector<cCalculator<double> *> aVEqCol ;

       for (size_t aKIm : aConfig)
       {
           aVS.push_back(mVSIm.at(aKIm));
           aVEqCol.push_back(mVSIm.at(aKIm)->GetEqColinearity());
       }

       const auto & aVals  = aPair.second;
       size_t aNbIm = aConfig.size();
       size_t aNbPts = aVals.mVPGround.size();


       //  parse all the multiple tie points of a given config
       for (size_t aKPts=0; aKPts<aNbPts ; aKPts++)
       {
           const cPt3dr & aPGr = aVals.mVPGround.at(aKPts);
           cSetIORSNL_SameTmp<tREAL8>  aStrSubst(aPGr.ToStdVector());
           //StdOut() << aPGr.x() << " " << aPGr.y() << " " << aPGr.z() << std::endl;

           size_t aNbEqAdded = 0;
           for (size_t aKIm=0 ; aKIm<aNbIm ; aKIm++)
           {
               const cPt2dr & aPIm =  aVals.mVPIm.at(aKPts*aNbIm+aKIm);
               cSensorImage* aSens = aVS.at(aKIm);

// StdOut() << "VISSSS " << aSens->IsVisibleOnImFrame(aPIm) << " " <<  aSens->IsVisible(aPGr) << "\n";

               if (aSens->IsVisibleOnImFrame(aPIm) && aSens->IsVisible(aPGr))
               {
                   cPt2dr aResidual  = aPIm-aSens->Ground2Image(aPGr);
                   tREAL8 aWeightImage =  aTieP_Weighter.SingleWOfResidual(aResidual);

                   // StdOut() << "RRRR " << aResidual << " W=" << aWeightImage << "\n";

                   cCalculator<double> * anEqColin =  aVEqCol.at(aKIm);

                   std::vector<double> aVObs = aPIm.ToStdVector();  // put Xim & Yim as observation
                   aSens->PushOwnObsColinearity(aVObs,aPGr);  // add eventual observation of sensor (as rot with central persp)

                   std::vector<int> aVIndGlob = {-1,-2,-3};  // index of unknown, begins with temporay
                   for (auto & anObj : aSens->GetAllUK())  // now put sensor unknown
                   {
                      anObj->PushIndexes(aVIndGlob);
                   }

                   tREAL8 aN2Res = Norm2(aResidual);
                   if (aWeightImage>0)
                   {
                       aWeigthedRes.Add(aWeightImage,aN2Res);
                       mSys->R_AddEq2Subst(aStrSubst,anEqColin,aVIndGlob,aVObs,aWeightImage);
                       aNbEqAdded++;
                   }
                   aUnWeightedRed.Add(1.0,aN2Res);
                   if (doStatRes)
                      aStatResidual.Add(aN2Res);
               }
           }

           // if at least 2 tie-point, we can add equation with schurr-complement
           if (aNbEqAdded>=2)
              mSys->R_AddObsWithTmpUK(aStrSubst,mCurLVMParam);  // finnaly add obs accummulated
       }

   }
   StdOut() << Color::argument << " *[TieP]: " << Color::end
            << Color::descr <<  "  Name=" << Color::end << aBA_TieP.mName
            << Color::descr << " Weighted="  << Color::end << aWeigthedRes.Average(-1)
            << Color::descr << " UnWeighted=" << Color::end << aUnWeightedRed.Average(-1)
            << Color::descr << " NbObs=" << Color::end  << aWeigthedRes.Nb()
     ;

   /// MPD :done quick & dirty,  amelioration : user fix %% + put in MMVII lib the format
   if (aStatResidual.NbMeasures()!=0)
   {
       StdOut()  << Color::descr << " ;; ResAt[x%%]=YY, {" << Color::end;
       for (const auto & aPerMil : mTiePShowPerMil)
       {
            tREAL8 aR = aStatResidual.ErrAtProp(aPerMil/1000.0) ;

            StdOut()  << Color::descr<<  aPerMil << ":"  << Color::end << ResidualToStr(aR)  << "  ";
       }
       StdOut() << Color::descr << "}" << Color::end;
   }

   StdOut() << std::endl;
}

void cMMVII_BundleAdj::OneItere_TieP()
{
    for (const auto & aParamTieP  : mVTieP)
        OneItere_TieP(*aParamTieP);
}


};



