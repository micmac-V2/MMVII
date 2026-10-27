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

         tREAL8 ScorePairs() const;
         tREAL8 ScoreTriplet() const;


         void AddPair(const std::string& aIm1,const std::string& aIm2);
         cPerspCamIntrCalib  * CalibOfFocalePixel(tREAL8 aFocalePix);
         void TestOneFoc35(tREAL8 aFocEq35);

         tRes1Pair  EstimatePose2IM(const std::string& aIm1,const std::string& aIm2);

        cPhotogrammetricProject   mPhProj;
        std::string               mPrefixTmp;


        // =========== Mandatory args ================
        std::string               mPatIm;

        // =========== Optionnal args ================

        eStepOriRelPipeline       mStepEnd;
        eProjPC                   mTypeProj;
        cPt2dr                    mIntervFocal35mm;
        tREAL8                    mRatioFocSucc;

        cPt2di                    mSzIm;
        std::string               mNameCalib;

        cPerspCamIntrCalib        * mCalib;
        cPt3di                    mDegree;




        cTimerSegm *              mTimeSegm ;
        bool                      mShow;
};

cAppli_OriRelEstimFocale::cAppli_OriRelEstimFocale(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
    cMMVII_Appli    (aVArgs,aSpec),
    mPhProj         (*this),
    mStepEnd        (eStepOriRelPipeline::eOriTriplets),
    mTypeProj       (eProjPC::eStenope),
    mRatioFocSucc   (1.1),
    mDegree         (3,1,1),
    mShow           (false)
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
            << AOpt2007(mStepEnd,NP_StepEndEstimFocal,"Step at which it eval is made",{eTA2007::HDV})


      ;

    return anArgOpt;
}

std::vector<cOneHelpSampleCmp>  cAppli_OriRelEstimFocale::Samples() const
{

    return {};
}

// MMVII OriPoseEstimRelGlob .*JPG Rel CalibInit InTieP=V1 VirTP=5Pt OutOri=GlobInit

tREAL8 cAppli_OriRelEstimFocale::ScorePairs() const
{
   // std::string OriRel_NamePairsOfAllImages(bool isIn, std::string aPost="" ) const;

    cWeightAv<tREAL8,tREAL8> aRes;

    std::vector<const tNamePair *>  aVecP;
    tNameRel  aSetOfPair;
    ReadFromFile(aSetOfPair, mPhProj.OriRel_NamePairsOfAllImages(false));
    aSetOfPair.PutInVect(aVecP,false);

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
    cWeightAv<tREAL8,tREAL8> aRes;
    tNameSet aSetN1;
    ReadFromFile(aSetN1,mPhProj.OriRel_NameAllImages(false));
    std::vector<std::string > aVecName = ToVect(aSetN1);


    std::vector<cDataSolOriTriplet>  aV3 = mPhProj.ReadAllTriplets(ToVect(aSetN1),false,false);

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
    tREAL8 aFocPix = aFocEq35 * Norm2(mSzIm) / Norm2(cPt2dr(24,36));
    cPerspCamIntrCalib  * aCalib = CalibOfFocalePixel(aFocPix );

    //std::string mPrefix = TheSpec_OriEstimFocale.Name();
    std::string aDir =     cDirsPhProj::StdDirOfType(eTA2007::Orient,mPrefixTmp);
    CreateDirectories(aDir);
    aCalib->ToFile(aDir+mNameCalib+"."+TaggedNameDefSerial());


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

    tREAL8 aScore2 = ScorePairs();
     tREAL8 aScore3 =  ((int)mStepEnd >= (int) eStepOriRelPipeline::eOriTriplets) ? ScoreTriplet() : 0.0;

    StdOut() << "Test F35 : " << aFocEq35
             << " => " << aFocPix
             << " SC2=" << aScore2
             << " SC3=" << aScore3
             << "\n";

 /*   if ((int)mStepEnd >= (int) eStepOriRelPipeline::eOriTriplets)
    {
           ScoreTriplet();
    }*/

   //  getchar();
    delete aCalib;

}





int cAppli_OriRelEstimFocale::Exe()
{
    mTimeSegm = mShow ? new cTimerSegm(this) : nullptr ;
    mPhProj.FinishInit();


    mPrefixTmp   = mPhProj.DPOriRel().DirOut();


    // ----- Estimate the interval of focal if not specified by user ------------------
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

    // ----------------  Compute size of camera/image, Name of calib and check = for all images ----------
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


    tREAL8 aLogF0 = std::log(mIntervFocal35mm.x());
    tREAL8 aLogF1 = std::log(mIntervFocal35mm.y());
    int aNbTest = round_up ((aLogF1-aLogF0)/ std::log(mRatioFocSucc));
    tREAL8 aLogStep = (aLogF1-aLogF0) / aNbTest;

    for (int aKTest=0 ; aKTest<=aNbTest ; aKTest++)
    {
        tREAL8 aLogF = aLogF0 + aKTest*aLogStep;
        tREAL8 aFoc35 = std::exp(aLogF);
        TestOneFoc35(aFoc35);
    }

    delete mTimeSegm;
    return EXIT_SUCCESS;
}


/* ====================================================== */
/*               OriPoseEstimRel2Im                       */
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




