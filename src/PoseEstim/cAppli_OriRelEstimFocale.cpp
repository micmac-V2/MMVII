#include "MMVII_TplHeap.h"

#include "MMVII_PoseRel.h"
#include "MMVII_Tpl_Images.h"
#include "MMVII_HeuristikOpt.h"
#include "MMVII_DeclareAllCmd.h"
#include "MMVII_Random.h"

namespace MMVII
{

const std::string  NP_CalledByEstimFocal="CalledByEstimFocal";
const std::string  NP_StepEndEstimFocal = "StepEnd";

/*
class cPairImEstimFocale
{
    public :
        cPairImEstimFocale(const std::string & aIm1,const std::string & aIm2);
        std::string mIm1;
        std::string mIm2;

        cSetHomogCpleIm           mCpleHFull;
        cSetHomogCpleIm           mCpleBig;
        cSetHomogCpleIm           mCpleAvg;
        cSetHomogCpleIm           mCpleSmall;
};

cPairImEstimFocale::cPairImEstimFocale(const std::string & aIm1,const std::string & aIm2) :
    mIm1 (aIm1),
    mIm2 (aIm2)
{
    OrderMinMax(mIm1,mIm2);
}
*/

   /* ********************************************************** */
   /*                                                            */
   /*                 cAppli_OriRel2Im                           */
   /*                                                            */
   /* ********************************************************** */

class cAppli_OriRelEstimFocale : public cMMVII_Appli
{
     public :

        typedef std::pair<int,cCdtFinalPoseRel2Im> tRes1Pair;
        static constexpr int RESULT_NO_POSE = 2007;

        typedef cIsometry3D<tREAL8>  tPose;

        cAppli_OriRelEstimFocale(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec);
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;
         std::vector<cOneHelpSampleCmp>  Samples() const override;

     private :

         /// return the score by reading the residual of pair orientation
         tREAL8 ScorePairs() const;

         /// return the score by reading the residual of triplet orientation
         tREAL8 ScoreTriplet() const;

          /// Compute the basic calibration w/o distorstion
         cPerspCamIntrCalib  * CalibOfFocalePixel(tREAL8 aFocalePix);

         /// make the test for a given focal, and register results
         void     TestOneFoc35(tREAL8 aFocEq35);

         /// Compute the focal lenght of given index
         tREAL8   FocOfLogInd(const tREAL8 anIndex) const;


        cPhotogrammetricProject   mPhProj;
        std::string               mPrefixTmp;  ///< Prefix used for temporary results (calib & OriRel)


        // =========== Mandatory args ================
        std::string               mPatIm;     ///< Images used

        // =========== Optionnal args ================

        eStepOriRelPipeline       mStepEnd;           ///< At which step of compute (pair, triplet ...) do we stop ?
        eProjPC                   mTypeProj;          ///< Type of camera projection
        cPt2dr                    mIntervFocal35mm;   ///<  Min-max focal to test
        tREAL8                    mRatioFocSucc;      ///< Ratio between to focal teste

        // ============ Internal parameter
        cPt2di                    mSzIm;      ///< size of images computed
        std::string               mNameCalib; ///< name of calibration
        cPt3di                    mDegree;    ///< Fix to 3,1,1 now

        int                      mKCurTest;  ///<  Num or current test
        cWhichMin<int,tREAL8>    mMinRes;    ///< structure to extract num test minimizing score
        std::vector<cPt2dr>      mVFocRes;   ///< Pair Foc/Residual :

        tREAL8                    mLogF0;       ///< log of min focal
        tREAL8                    mLogF1;       ///< log of max focal
        tREAL8                    mLogStep;     ///< log of each ratio between succ focal
        int                       mNbFocTested; ///< number of focal to  teste

};

cAppli_OriRelEstimFocale::cAppli_OriRelEstimFocale(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
    cMMVII_Appli    (aVArgs,aSpec),
    mPhProj         (*this),
    mStepEnd        (eStepOriRelPipeline::eOriTriplets),
    mTypeProj       (eProjPC::eStenope),
    mRatioFocSucc   (1.3),
    mDegree         (3,1,1),
    mKCurTest       (0)
{
}


cCollecSpecArg2007 & cAppli_OriRelEstimFocale::ArgObl(cCollecSpecArg2007 & anArgObl)
{

     anArgObl
             << Arg2007(mPatIm,"Pattern of images",{{eTA2007::MPatFile,"0"},{eTA2007::FileDirProj}})
             << mPhProj.DPOriRel().ArgDirOutMand()
     ;


     return anArgObl;
}

//cCollecSpecArg2007 & operator << (const cHeaderSectionArg& aComment);



cCollecSpecArg2007 & cAppli_OriRelEstimFocale::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    anArgOpt
            << cHeaderSectionArg("Input for tie points (at least one must be provided)")
            << mPhProj.AddAllDpTieP(anArgOpt)

            << cHeaderSectionArg("Parameter for calibration")
            << AOpt2007(mIntervFocal35mm,"IntFoc35","Interval focal in equiv 35mm, def value depend of model")
            << AOpt2007(mTypeProj,"Proj","Type of projection",{eTA2007::HDV})
            << AOpt2007(mRatioFocSucc,"RatioFocal","Ratio between sucessively tested focals",{eTA2007::HDV})
            << AOpt2007(mNbFocTested,"NbStep","Fix number of step, rather than ratio",{eTA2007::HDV})

            //  Supress 4 now acces to this option : Pairs in no good, Global is not implemented
            //  << AOpt2007(mStepEnd,NP_StepEndEstimFocal,"Step of process at which it eval is made",{eTA2007::HDV})
      ;

    return anArgOpt;
}

std::vector<cOneHelpSampleCmp>  cAppli_OriRelEstimFocale::Samples() const
{

    return {
               cOneHelpSampleCmp::Header("Examples with Ramses data set, in two steps"),
               {
                    "MMVII OriPoseEstimFocale \"IMG_034[036].JPG\" TMP InTieP=V1",

               },
               {
                   "MMVII OriPoseEstimFocale \"IMG_034[036].JPG\" TMP InTieP=V1  IntFoc35=[40,60] RatioFocal=1.1"
               },

               cOneHelpSampleCmp::Header("Examples with Niche data set"),
               {
                   "MMVII OriPoseEstimFocale \"IMG_569[456].JPG\" TMP InTieP=V1  Proj=FE_EquiDist"
               }
    };
}


tREAL8 cAppli_OriRelEstimFocale::ScorePairs() const
{
    // extract pairs computed
    std::vector<const tNamePair *>  aVecP;
    tNameRel  aSetOfPair;
    ReadFromFile(aSetOfPair, mPhProj.OriRel_NamePairsOfAllImages(false));
    aSetOfPair.PutInVect(aVecP,false);

    // compute average of score
    cWeightAv<tREAL8,tREAL8> aRes;
    for (const auto  & aPair : aVecP )
    {
        std::string aNameFile = mPhProj.OriRel_NameOriPair2Images(aPair->V1(),aPair->V2(),false);
        cCdtFinalPoseRel2Im aCdt;
        ReadFromFile(aCdt,aNameFile);
        aRes.Add(1.0,aCdt.mVCdt.at(0).mScore);
    }

    return aRes.Average();
}

tREAL8 cAppli_OriRelEstimFocale::ScoreTriplet() const
{
    // extract triplet computed
    tNameSet aSetN1;
    ReadFromFile(aSetN1,mPhProj.OriRel_NameAllImages(false));
    std::vector<std::string > aVecName = ToVect(aSetN1);
    std::vector<cDataSolOriTriplet>  aV3 = mPhProj.ReadAllTriplets(ToVect(aSetN1),false,false);

    // compute average of scores
    cWeightAv<tREAL8,tREAL8> aRes;
    for (const auto & a3 : aV3)
    {
        aRes.Add(1.0,a3.mScore);
    }

    return aRes.Average();
}


cPerspCamIntrCalib  * cAppli_OriRelEstimFocale::CalibOfFocalePixel(tREAL8 aFocalePix)
{
    return   cPerspCamIntrCalib::SimpleCalib
             (
                  mNameCalib,
                  mTypeProj,
                  mSzIm,
                  cPt3dr(mSzIm.x()/2.0,mSzIm.y()/2.0,aFocalePix),
                  mDegree
            );
}


void cAppli_OriRelEstimFocale::TestOneFoc35(tREAL8 aFocEq35)
{
    // ----------- Generate file with calibration for use by command
    tREAL8 aFocPix = F35toFPix(aFocEq35,mSzIm) ; // aFocEq35 * Norm2(mSzIm) / Norm2(cPt2dr(24,36));
    cPerspCamIntrCalib  * aCalib = CalibOfFocalePixel(aFocPix );
    std::string aDir =     cDirsPhProj::StdDirOfType(eTA2007::Orient,mPrefixTmp);
    CreateDirectories(aDir);
    aCalib->ToFile(aDir+mNameCalib+"."+TaggedNameDefSerial());


    //   ------------ Generate and execut the command -----------------
    {
        cParamCallSys aComPipeline(TheSpec_OriPoseEstimRel,true);
        aComPipeline.AddArgs(mPatIm);      // images
        aComPipeline.AddArgs(mPrefixTmp);  // Ori Rel
        aComPipeline.AddArgs(mPrefixTmp);  // Ori Calib

        aComPipeline.AddMMVIIAllTiepInInit(mPhProj);
        aComPipeline.AddMMVIIArgsOpt(NP_CalledByEstimFocal,"true");
        aComPipeline.AddMMVIIArgsOpt(NP_StepEndEstimFocal,E2Str(mStepEnd));

        aComPipeline.ExecuteWithoutMsg();
    }

    // ----------------- Extract the score, print, and update structures
    tREAL8 aScore2 = ScorePairs();
    tREAL8 aScore3 =  ((int)mStepEnd >= (int) eStepOriRelPipeline::eOriTriplets) ? ScoreTriplet() : 0.0;

    StdOut() << Color::title << " * OneTest, "
             << Color::descr << " F35= "    << Color::end  <<   ClaudeToFixDigit(aFocEq35,5)
             << Color::descr << "  FPix= "  << Color::end  <<   ClaudeToFixDigit(aFocPix,5)
             << Color::descr << "  SC2="    << Color::end  <<   ClaudeToFixDigit(aScore2,5)
             << Color::descr << "  SC3="    << Color::end  <<   ClaudeToFixDigit(aScore3,5)
             << "\n";

    tREAL8 aScore = aScore3;
    mMinRes.Add(mKCurTest,aScore);
    mVFocRes.push_back(cPt2dr(aFocEq35,aScore));

    // ------------ final
    delete aCalib;
    mKCurTest++;
}


tREAL8   cAppli_OriRelEstimFocale::FocOfLogInd(const tREAL8 aKTest) const
{
    return  std::exp(mLogF0 + aKTest*mLogStep);
}



int cAppli_OriRelEstimFocale::Exe()
{
    mPhProj.FinishInit();

    mPhProj.AssertAtLeastOneTiePModeIsInit();

    mPrefixTmp   = mPhProj.DPOriRel().DirOut();


    // -----  [1]  Estimate the interval of focal if not specified by user ------------------
    if (!IsInit(&mIntervFocal35mm))
    {
        switch (mTypeProj) {

            case eProjPC::eStenope :
               mIntervFocal35mm = cPt2dr(15,800);
            break;

            case eProjPC::eFE_EquiDist :
            case eProjPC::eFE_EquiSolid :
            case eProjPC::eStereroGraphik :
            case eProjPC::eOrthoGraphik :
                mIntervFocal35mm = cPt2dr(6,25);
            break;

            case eProjPC::eEquiRect :
                 mIntervFocal35mm = cPt2dr(3,15);
            break;

            default :
                 MMVII_UnclasseUsEr("Unhandled type of projection");
            break;
        }

        StdOut() << " Default interval of focal is : " <<  mIntervFocal35mm <<"\n";
    }


    std::vector<std::string> aVNamesIm =  VectMainSet(0);


    // ------------- [2]   Compute size of camera/image, Name of calib and check = for all images ----------
    for (size_t aK1 =0 ; aK1<aVNamesIm.size(); aK1++)
    {
        std::string aNameIm1 = aVNamesIm.at(aK1);


        std::string aNameCalib = mPhProj.StdNameCalibOfImage(aNameIm1);
       // StdOut() << " GGGG=" << aNameCalib << "\n";

        cDataFileIm2D aDFI = cDataFileIm2D::Create(aNameIm1,eForceGray::No);
        cPt2di aSzIm = aDFI.Sz();
        if (aK1==0)
        {
            mSzIm = aSzIm;
            mNameCalib = aNameCalib;
        }
        if (mSzIm!= aSzIm)
        {
            StdOut() << Color::error <<  " Sizes got " << Color::end << mSzIm << " "  << aSzIm << "\n";
            MMVII_UnclasseUsEr
            (
                  "Different image size for : " +  aVNamesIm.at(0) + " " + aNameIm1
            );
        }
        if (aNameCalib!= mNameCalib)
        {
            StdOut() << Color::error <<  " Names got " << Color::end << mNameCalib << " "  << aNameCalib << "\n";

            MMVII_UnclasseUsEr
            (
                  "Different calib name for : " +  aVNamesIm.at(0) + " " + aNameIm1
            );
        }
    }

    // ------------- [3]  Test the differents focal  ----------

    mLogF0 = std::log(mIntervFocal35mm.x());
    mLogF1 = std::log(mIntervFocal35mm.y());
    if (!IsInit(&mNbFocTested))
        mNbFocTested = round_up ((mLogF1-mLogF0)/ std::log(mRatioFocSucc));
    mLogStep = (mLogF1-mLogF0) / mNbFocTested;

    for (int aKTest=0 ; aKTest<=mNbFocTested ; aKTest++)
    {
       TestOneFoc35(FocOfLogInd(aKTest));
    }

    //------------ [4]  Print the result -------------------------------

             // [4.1]  print "basic" min

    int aKMin = mMinRes.IndexExtre();
    cPt2dr aFRoc = mVFocRes.at(aKMin);

    StdOut() << "\n" << Color::title  <<"============ Final  estimation =========== \n\n";
    StdOut() <<  "   "
             << Color::descr << "F35= "    << Color::end  <<   ClaudeToFixDigit(aFRoc.x(),5)
              << Color::descr << " Res= "    << Color::end  <<   ClaudeToFixDigit(aFRoc.y(),5)
             << "\n";

             // [4.2]  Extract and print parabolic refinement

    if ((aKMin>=1) && (aKMin<((int)mVFocRes.size()-1)))
    {

          tREAL8 aS0 =  mVFocRes.at(aKMin-1).y();
          tREAL8 aS1 =  mVFocRes.at(aKMin  ).y();
          tREAL8 aS2 =  mVFocRes.at(aKMin+1).y();

           tREAL8 aSIE = StableInterpoleExtr(aS0,aS1,aS2);

           tREAL8 aF35Refine =  FocOfLogInd(aKMin+aSIE);
           tREAL8 aFPixRefine = F35toFPix(aF35Refine,mSzIm);
           StdOut() << Color::success << "     --- F35-Parabolic-Refine ---- " << "\n";
           StdOut() <<  "   "
                    << Color::descr << "F35= "     << Color::end  <<   ClaudeToFixDigit(aF35Refine,5)
                    << Color::descr << " FPix= "    << Color::end  <<  ClaudeToFixDigit(aFPixRefine,5);

           cMetaDataImage  aMTD = mPhProj.GetMetaData(aVNamesIm.at(0)) ;
           const cElemCamDataBase *  anEDB = mPhProj.GetCamFromNameCam(aMTD.CameraName(),true);
           if (anEDB)
           {
               tREAL8 aFocReal = aF35Refine * Norm2(anEDB->mSzSensor_Mm) / Norm2(cPt2dr(36,24)) ;
               StdOut() << Color::descr << " Freal= "     << Color::end  <<   ClaudeToFixDigit(aFocReal,5);
           }
           StdOut ()      << "\n";

    }
    else
    {
        StdOut() << "\n" << Color::warning  <<" Min val for extremum : no parabolic refinement \n";
    }


    StdOut() << "\n" << Color::title  <<"========================================== \n" << Color::end;

    return EXIT_SUCCESS;
}


/* ====================================================== */
/*               OriPoseEstimFocale                       */
/* ====================================================== */

tMMVII_UnikPApli Alloc_OriRelEstimFocale(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppli_OriRelEstimFocale(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_OriEstimFocale
(
     "OriPoseEstimFocale",
      Alloc_OriRelEstimFocale,
      "Rough estimation of focale length",
      {eApF::Ori},
      {eApDT::TieP},
      {eApDT::Orient},
      __FILE__
);


}; // MMVII




