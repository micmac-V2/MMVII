#include "MMVII_InstrumentalBlock.h"
#include "cMMVII_Appli.h"
#include "MMVII_2Include_Serial_Tpl.h"
#include "MMVII_HeuristikOpt.h"
#include "MMVII_Tpl_OptimManifold.h"


/**
  \file cInstrumentalBloc.cpp


  \brief This file contains the core implemantation of Block of rigid instrument
 
*/

namespace MMVII
{

static const std::string N_Indep     = "D2_Indep";
static const std::string N_2Orthog   = "D3_2Orthog";
static const std::string N_VertI     = "D4_VertIndep";
static const std::string N_Vert2O    = "D5_Vert2Orthog";



/* *************************************************************** */
/*                                                                 */
/*               cAppli_BlockInstrInitClino                        */
/*                                                                 */
/* *************************************************************** */

class cRes_ClinoInit
{
   public :
       cRes_ClinoInit();

       cRes_ClinoInit(size_t aNb,size_t aDOF,bool forceOrthog,bool isVertFree,const std::string & aMsg);

       void SetKthClino(int aK,tREAL8 aScore,const cPt3dr &);
       void AddVertical(const cPt3dr &);


       size_t                    mDOF;
       bool                      mForceOrthog;
       bool                      mVertFree;
       std::string               mMsg;
       cWeightAv<tREAL8,tREAL8>  mAvg;
       std::vector<tREAL8>       mVScore;
       std::vector<cPt3dr>       mVDir;
       std::vector<cPt3dr>       mVVert;
        cVarPts<3>               mVarVert;
};

cRes_ClinoInit::cRes_ClinoInit()
{
}

cRes_ClinoInit::cRes_ClinoInit(size_t aNb,size_t aDOF,bool forceOrthog,bool isVertFree,const std::string & aMsg) :
   mDOF          (aDOF),
   mForceOrthog  (forceOrthog),
   mVertFree     (isVertFree),
   mMsg          (aMsg),
   mVScore       (aNb,-1.0),
   mVDir         (aNb,cPt3dr(0,0,0))
{
}

void cRes_ClinoInit::AddVertical(const cPt3dr & aVert)
{
    mVVert.push_back(aVert);
    mVarVert.Add(aVert);
}


void cRes_ClinoInit::SetKthClino(int aK,tREAL8 aScore,const cPt3dr & aDir)
{
    MMVII_INTERNAL_ASSERT_always(mVScore.at(aK)<0,"Multiple Init in cRes_ClinoInit");

    mAvg.Add(1.0,aScore);
    mVScore.at(aK) = aScore;
    mVDir.at(aK)   = aDir;
}

/**
 * @brief The cAppli_BlockInstrInitClino class
 *
 * Class for initializing the calibration of clino from a set images oriented with clino measurement
 */




class cAppli_BlockInstrInitClino : public   cMMVII_Appli,
                                   public   cOptDiscScorer<cPt3dr>,           // Compute 1 Clino
                                   public   cOptDiscScorer<tRotR>,             // Compute 2 Clino Orthog
                                   public   cOptDiscScorer<tPairPt3dr>,  // Compute Vert + 1 Clino
                                   public   cOptDiscScorer<tPoseAsPair>        // Compute Vert + 2 Clino Orthog
{
     public :
        // ===============   method to make it a full cMMVII_Appli ========
        cAppli_BlockInstrInitClino(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli &);
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;
        int Exe() override;
        // std::vector<std::string>  Samples() const ;
     private :

        /// Angle in the unit selected for print (DmGon in general)
        tREAL8 Ang4Show(tREAL8) const;

        // =============  Methods for individual clino ====================

        tREAL8 Score(const cPt3dr&) const override;

           /// Calibrate a single clino independantly of others
        void Process1ClinoIndep(int aK);
           /// Score of calibration "aDir"  for clino "aKClino"
        tREAL8 Score1Clino(const cPt3dr& aDir,size_t aKClino) const;
           /** Interface to  ScoreDirClino for cDataMapping  R2 -> R ,
            once it has been roughly initialized compute in tangent plane */



        // =============  Methods for orthohonal clino ====================

        /*  If we have 2 clino, forced to be orthogonal, we modelize them with
         *  a rotation R such axe I and J of the rotation are the axes of clino 1 and 2
         */
        tREAL8 Score(const tRotR&) const override;

        /// Calibrate a pair of clino forced to be orthog
        void Process2OrthogClino(const cPt2di &);

        /// Score of pair of orthog clino coded as a rotation
          tREAL8 Score2OthogClino(const tRotR& aR) const;
        /// Score of an WPK, with signature for cDataMapping
       ///  cPt1dr  Value(const cPt3dr&) const override ;

          void AddResult(const std::string&,int aDOF,bool forceOrthog,bool isVertFree);
          bool  ResultIsInit(const std::string&) const;
          cRes_ClinoInit & Result(const std::string&);

          // =============  Methods for 1 clino + Vertical ====================

          tREAL8 Score(const tPairPt3dr&) const override;
          void ProcessVerticalAnd1Clino(int  aKClino);


        // =============  Methods for 2 orthogonal clino + Vertical ====================

          tREAL8 Score(const tPoseAsPair&) const override;
          void ProcessVerticalAnd2ClinoOrthog(const cPt2di & aK1K2);


        cPhotogrammetricProject   mPhProj;


        // -------------------------- Mandatory parameters ---------------------------
        std::string               mSpecImIn;  ///< images
        std::vector<int>          mDOFTested;     ///< 5-> Vert and 2, 4 -> Vert & Indep, 3-> Orthog, 2 ->Indep

        // ------------- Optionnal parameters ----------------
        std::string               mNameBloc;     //< name of the bloc inside the
        int                       mNbSamp1Cl;     //< Number of sample of sphere for 1 clino (in each face of cube)
        cPt2di                    mNbSVert;      //< Nb Sample Vert&Rot, in case we guess vertical
        std::vector<cPt2di>       mPairOrthog;   //< vector of orthog if we change (more for test)
        eTyUnitAngle              mUnitShow;     //< Unity for print result
        bool                      mReSetSigma;   //< do we use sigma

        //  ----------
        cIrbComp_Block *          mBlock;        //< Bloc : Calib + Data of time stamp
        cIrbCal_Block *           mCalBlock;     //< Calibration of above block
        cIrbCal_ClinoSet *        mSetClinos ;   //< Set of clino of above calib
        int                       mNbClino;      //< Number of clinos



        int                       mKCurClino;     //< Index of cur clino, fix by DoOneClino
        int                       mK1CurC;        //< Case we have 2 clino, Index of Clin1
        int                       mK2CurC;        //<                      , Index of Clin2

        std::vector<int>          mNumPoseInstr;  //< Num cams used for pose estimation of instrument
        std::map<std::string,cRes_ClinoInit>  mMapRes;
        cRes_ClinoInit *                      mLastResult;

        /*
        cWeightAv<tREAL8,tREAL8>  mAvgClinIndep;
        std::vector<tREAL8>       mScoreClinoIndep;
        std::vector<cPt3dr>       mDirClinoIndep;

        cWeightAv<tREAL8,tREAL8>  mAvgClinGlob;
        std::vector<tREAL8>       mScoreClinoGlob;
        std::vector<cPt3dr>       mDirClinoGlob;
        */

        cPt3dr                    mVertApriori;
/*
        bool                      mDoVertAnd2C;
        bool                      mDoVertAnd1C;
        bool                      mDo2C;
        bool                      mDo1C;
        */
};



cAppli_BlockInstrInitClino::cAppli_BlockInstrInitClino
(
        const std::vector<std::string> &  aVArgs,
        const cSpecMMVII_Appli & aSpec
) :
    cMMVII_Appli (aVArgs,aSpec),
    mPhProj        (*this),
    mNameBloc      (cIrbCal_Block::theDefaultName),
    mNbSamp1Cl     (100),
    mNbSVert       (7,7),
    mUnitShow      (eTyUnitAngle::eUA_DMgon),
    mReSetSigma    (true),
    mBlock         (nullptr),
    mCalBlock      (nullptr),
    mSetClinos     (nullptr),
    mNbClino       (-1),
    mKCurClino     (-1),
    mLastResult    (nullptr),
    mVertApriori   (0,0,1)
{
}


cCollecSpecArg2007 & cAppli_BlockInstrInitClino::ArgObl(cCollecSpecArg2007 & anArgObl)
{
     return anArgObl
             <<  Arg2007(mSpecImIn,"Pattern/file for images", {{eTA2007::MPatFile,"0"},{eTA2007::FileDirProj}}  )
              << Arg2007(mDOFTested,"Degree of freedom tested in 2,3,4,5 ",{{eTA2007::ISizeV,"[1,4]"}})
             <<  mPhProj.DPBlockInstr().ArgDirInMand()
             <<  mPhProj.DPOrient().ArgDirInMand()
             <<  mPhProj.DPMeasuresClino().ArgDirInMand()
             <<  mPhProj.DPBlockInstr().ArgDirOutMand()
     ;
}

cCollecSpecArg2007 & cAppli_BlockInstrInitClino::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
        return anArgOpt
            << AOpt2007(mNameBloc,"NameBloc","Name of bloc to calib ",{{eTA2007::HDV}})
            << AOpt2007(mNbSamp1Cl,"NbSS","Number of sample on the sphere ",{{eTA2007::HDV}})
            << AOpt2007(mNumPoseInstr,"NPI","Num of cams used  for estimate pose of intsrument")
            << AOpt2007(mPairOrthog,"PairOrthog","Pair of orthogonal camera (if reset)")
            << AOpt2007(mUnitShow,"USA","Unity Show Angles",{AC_ListVal<eTyUnitAngle>(),{eTA2007::HDV}})
            << AOpt2007(mReSetSigma,"ResetSigma","Do we use sigma",{{eTA2007::HDV}})
            << AOpt2007(mVertApriori,"Vert0","A priori vertical => we compte both clino & vert")
            << AOpt2007(mNbSVert,"NbSV","Number of sample, when both Orthog & Vert, x=>quat, y=>Spher",{eTA2007::HDV})
    ;
}


tREAL8 cAppli_BlockInstrInitClino::Ang4Show(tREAL8 anAng) const
{

   return AngleFromRad(anAng,mUnitShow);
}

void cAppli_BlockInstrInitClino::AddResult(const std::string & aName,int aDOF,bool forceOrthog,bool isVertFree)
{
    mMapRes[aName] = cRes_ClinoInit(mNbClino,aDOF,forceOrthog,isVertFree,aName);
    mLastResult = &  mMapRes[aName];
}

bool  cAppli_BlockInstrInitClino::ResultIsInit(const std::string& aName) const
{
    return MapBoolFind(mMapRes,aName);
}

cRes_ClinoInit & cAppli_BlockInstrInitClino::Result(const std::string& aName)
{
    MMVII_INTERNAL_ASSERT_medium(ResultIsInit(aName)," cAppli_BlockInstrInitClino::Result");
    return mMapRes[aName] ;
}


    //=====================================================================
    // ========================  5D optimisation ==========================
    //=====================================================================


tREAL8  cAppli_BlockInstrInitClino::Score(const tPoseAsPair& aPair) const
{
    mBlock->SetVerticalCste(VUnit(aPair.first));
    return Score2OthogClino(aPair.second);

}

void cAppli_BlockInstrInitClino::ProcessVerticalAnd2ClinoOrthog(const cPt2di & aK1K2)
{
    mK1CurC = aK1K2.x();
    mK2CurC = aK1K2.y();

    // Let  K1,K2 be two direction of clino and V a vertical, we have the equivalence
    //                (K1,K2,V)  ~ (-K1,-K2,-V)  =>  true : ambiguity above
    tPoseCart_OptimDisc aPoseDesc = Pose_OptimDisc(mNbSVert.y(),true,mNbSVert.x());

    // object for optimizing
    cTplOptDisc_OnManifold<tPoseCart_OptimDisc> anOptim(aPoseDesc,*this);

    anOptim.ComputeSol(0.5,1e-7,5);

    const cWhichMin<tPoseAsPair,tREAL8> & aSol = anOptim.GetSol();
    auto [aVert,aRot] =  aSol.IndexExtre();

    Result(N_Vert2O).SetKthClino(mK1CurC,aSol.ValExtre(),aRot.AxeI());
    Result(N_Vert2O).SetKthClino(mK2CurC,aSol.ValExtre(),aRot.AxeJ());

    Result(N_Vert2O).AddVertical(aVert);

    mBlock->SetVerticalCste(cPt3dr(0,0,1));
}

    //=====================================================================
    // ========================  4D optimisation ==========================
    //=====================================================================


tREAL8  cAppli_BlockInstrInitClino::Score(const tPairPt3dr& aPair) const
{
    mBlock->SetVerticalCste(VUnit(aPair.first));
    return Score1Clino(aPair.second,mKCurClino);
}

void cAppli_BlockInstrInitClino::ProcessVerticalAnd1Clino(int  aKClino)
{
    mKCurClino = aKClino;

    cTplOptDisc_OnManifold<t2Sph3_OptimDisc> anOptim(TwoSph3_OptimDisc(10,true,10,false),*this);
    anOptim.ComputeSol(0.5,1e-8,5);
    const cWhichMin<tPairPt3dr,tREAL8> & aSol = anOptim.GetSol();
    auto [aVert,aDirC] = aSol.IndexExtre();

    Result(N_VertI).SetKthClino(mKCurClino,aSol.ValExtre(),aDirC);
    Result(N_VertI).AddVertical(aVert);

    mBlock->SetVerticalCste(cPt3dr(0,0,1));

//    StdOut() << " ProcessVerticalAnd1Clino " << Ang4Show(aSol.ValExtre()) << " V=" << aSol.IndexExtre().first << "\n";
}

     //=====================================================================
     // ========================  3D optimisation ==========================
     //=====================================================================

tREAL8 cAppli_BlockInstrInitClino::Score2OthogClino(const tRotR& aR) const
{
        tREAL8 aRes1 = Score1Clino(aR.AxeI(),mK1CurC);
        tREAL8 aRes2 = Score1Clino(aR.AxeJ(),mK2CurC);

        return (aRes1+aRes2)/2.0;
}

tREAL8 cAppli_BlockInstrInitClino::Score(const tRotR& aR) const
{
        return Score2OthogClino(aR);
}


void cAppli_BlockInstrInitClino::Process2OrthogClino(const cPt2di & aK1K2)
{
// StdOut() <<  "Process2OrthogClinoProcess2OrthogClino " << aK1K2 << "\n";

    mK1CurC = aK1K2.x();
    mK2CurC = aK1K2.y();

    cTplOptDisc_OnManifold<cRot_OptimDisc> anOptim(cRot_OptimDisc(15),*this);
    anOptim.ComputeSol(0.5,1e-8,5);
    const cWhichMin<tRotR,tREAL8> & aSol = anOptim.GetSol();

    Result(N_2Orthog).SetKthClino(mK1CurC,aSol.ValExtre(),aSol.IndexExtre().AxeI());
    Result(N_2Orthog).SetKthClino(mK2CurC,aSol.ValExtre(),aSol.IndexExtre().AxeJ());

   //   StdOut() << " RRRRRR " << Ang4Show(aSol.ValExtre()) << "\n";
}

            //=====================================================================
            // ========================  2D optimisation ==========================
            //=====================================================================


tREAL8 cAppli_BlockInstrInitClino::Score1Clino(const cPt3dr& aDir,size_t aKClino) const
{
        return mBlock->ScoreDirClino(aDir,aKClino);
}

tREAL8 cAppli_BlockInstrInitClino::Score(const cPt3dr& aDir) const
{
    return Score1Clino(aDir,mKCurClino);
}



void cAppli_BlockInstrInitClino::Process1ClinoIndep(int aK)
{
    mKCurClino = aK;


    cTplOptDisc_OnManifold<cSph3_OptimDisc> anOptim(cSph3_OptimDisc(50,false),*this);
    anOptim.ComputeSol(0.5,1e-8,5);
    const cWhichMin<cPt3dr,tREAL8> & aSol = anOptim.GetSol();


    Result(N_Indep).SetKthClino(mKCurClino,aSol.ValExtre(),aSol.IndexExtre());

}


            //=====================================================================
            // ===========================  ::Exe  ================================
            //=====================================================================

int cAppli_BlockInstrInitClino::Exe()
{
    mPhProj.FinishInit();


    // read an existing bloc from std folder
    mBlock = new cIrbComp_Block(mPhProj,mNameBloc);
    mCalBlock  = & mBlock->CalBlock();
    mSetClinos = & mCalBlock->SetClinos();
    mNbClino = mSetClinos->NbClino() ;

    // Initialize the structure for results

    for (const auto & aDOF : mDOFTested)
    {
        if      (aDOF==2)  AddResult(N_Indep,2,false,false);
        else if (aDOF==3)  AddResult(N_2Orthog,3,true,false);
        else if (aDOF==4)  AddResult(N_VertI,4,false,true);
        else if (aDOF==5)  AddResult(N_Vert2O,5,true,true);
        else
        {
            MMVII_UnclasseUsEr("Bad value in degree of freedom");
        }
    }

  //  mScoreClinoIndep.resize(mNbClino);
  //  mDirClinoIndep.resize(mNbClino);

    if (IsInit(&mNumPoseInstr))
    {
       mCalBlock->SetCams().SetNumPoseInstr(mNumPoseInstr);
    }


    //  add all the camera
    for (const auto & aNameIm :  VectMainSet(0))
    {
       mBlock->AddImagePose(aNameIm);
    }
    // Init vertical
    mBlock->SetVerticalCste(cPt3dr(0,0,1));



    mBlock->ComputePoseInstrument();
    mBlock->SetClinoValues();

    if (!IsInit(&mPairOrthog))
    {
        for (const auto & [aPair,aCstr] : mCalBlock->CstrOrthog())
        {
            int aK1 = mSetClinos->IndexClinoFromName(aPair.V1(),true);
            int aK2 = mSetClinos->IndexClinoFromName(aPair.V2(),true);

            if ((aK1>=0) && (aK2>=0))
                mPairOrthog.push_back(cPt2di(aK1,aK2));
        }
    }

    if (ResultIsInit(N_Vert2O))
    {
        for (const auto aPt: mPairOrthog)
        {
           ProcessVerticalAnd2ClinoOrthog(aPt);
        }
    }
    if (ResultIsInit(N_VertI))
    {
        for (int aKC=0 ; aKC<(int)mSetClinos->NbClino() ; aKC++)
             ProcessVerticalAnd1Clino(aKC);
    }
    if (ResultIsInit(N_2Orthog))
    {
        for (const auto aPt: mPairOrthog)
        {
            Process2OrthogClino(aPt);
        }
    }
    if (ResultIsInit(N_Indep))
    {
        for (int aKC=0 ; aKC<(int)mSetClinos->NbClino() ; aKC++)
            Process1ClinoIndep(aKC);
    }


    for (const auto & [aName,aRes] : mMapRes)
    {
        StdOut() << "=================== "  << aName
                << " AvgRes=" << Ang4Show(aRes.mAvg.Average())
                << " DMgon =====================\n";

        if (aRes.mVertFree)
        {
            StdOut() << "  * VERTICAL Avg=" << aRes.mVarVert.Avg() ;
            if ( aRes.mVVert.size() > 1)
                  StdOut()  << " StdDev=" << Ang4Show(aRes.mVarVert.UB_StdDev()) << " DMgon ";
             StdOut() << "\n";
        }

        if (! aRes.mForceOrthog)
        {
            StdOut() << "  * ORTHOGONALITY = ";
            for (const auto & aPt : mPairOrthog)
            {
                cPt3dr aD1 = aRes.mVDir.at(aPt.x());
                cPt3dr aD2 = aRes.mVDir.at(aPt.y());

                tREAL8 anAng = std::abs(AbsAngleTrnk( aD1,aD2)-M_PI/2.0) ;
                StdOut() << " " << aPt << Ang4Show(anAng) ;
                // StdOut() << " DDD " << aD1 << aD2 << Scal(aD1,aD2) ;
            }
            StdOut() << " DMgon \n";
        }

        StdOut() << "\n";
    }

    if (mLastResult!=nullptr)
    {
        StdOut() << " ------------ SAVE RES FOR : " << mLastResult->mMsg << " ------------\n";

        // ========================= print the result  (Reports 2 add) ===========================================
        for (int aKC=0 ; aKC<mNbClino ; aKC++)
        {
            tREAL8 aScore = mLastResult->mVScore.at(aKC);
            cPt3dr aDir = mLastResult->mVDir.at(aKC);

            std::string aName = mSetClinos->VNames().at(aKC);
            if (mReSetSigma)
            {
                auto & aDescInstr = mCalBlock->AddSigma_Indiv(aName,eTyInstr::eClino);
                //                                  WTr   WRot  STr   SRot
                aDescInstr.SetSigma(cIrb_SigmaInstr(0.0  ,1.0  ,0.0  ,aScore));
            }

            mSetClinos->KthClino(aKC).SetPNorm(aDir);

            StdOut() << " K=" << aKC << " N=" << aName
                     << " Score=" << Ang4Show(aScore)
                     << " Dir=" <<  aDir
                     << "\n";

        }
        mPhProj.SaveRigBoI(*mCalBlock);

        if (mLastResult->mVertFree)
        {
           mPhProj.DPOrient().SetDirOutInIfNotInit();

           mBlock->PtrOriSysco()->SetVertical(mLastResult->mVarVert.Avg());
           mPhProj.SaveCurSysCoOri( mBlock->PtrOriSysco());
        }
    }

    delete mBlock;
    return EXIT_SUCCESS;
}

    /* ==================================================== */
    /*                                                      */
    /*               MMVII                                  */
    /*                                                      */
    /* ==================================================== */


tMMVII_UnikPApli Alloc_BlockInstrInitClino(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppli_BlockInstrInitClino(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_BlockInstrInitClino
(
     "BlockInstrInitClino",
      Alloc_BlockInstrInitClino,
      "Init  camera poses inside a block of instrument",
      {eApF::BlockInstr,eApF::Ori,eApF::Clino},
      {eApDT::BlockInstr,eApDT::Ori,eApDT::Clino},
      {eApDT::BlockInstr},
      __FILE__
);


};

