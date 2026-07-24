#include "BundleAdjustment.h"
#include "MMVII_util_tpl.h"

#include "MMVII_Topo.h"

/**
   \file cAppliBundAdj.cpp

*/


namespace MMVII
{

//#define  TOTO(A,B) {StdOut() << A << B;}
//#define TOTO(A) {TOTO(A,"B")}

/* ************************************************************************ */
/*                                                                          */
/*                        cStdWeighterResidual                              */
/*                                                                          */
/* ************************************************************************ */

cStdWeighterResidual::cStdWeighterResidual(tREAL8 aSGlob,tREAL8 aSigAtt,tREAL8 aThr,tREAL8 aExp) :
   mWGlob       (1/Square(aSGlob)),
   mWithAtt     (aSigAtt>0),
   mSig2Att     (Square(aSigAtt)),
   mWithThr     (aThr>0),
   mSig2Thrs    (Square(aThr)),
   mExpS2       (aExp/2.0)
{
}

cStdWeighterResidual::cStdWeighterResidual(const std::vector<tREAL8> & aVect,int aK0) :
    cStdWeighterResidual
    (
        aVect.at(aK0),
        GetDef(aVect,aK0+1,-1.0),
        GetDef(aVect,aK0+2,-1.0),
        GetDef(aVect,aK0+3, 1.0)
    )
{
}

cStdWeighterResidual::cStdWeighterResidual() :
    cStdWeighterResidual({1.0},0)
{
}


tREAL8  cStdWeighterResidual::SingleWOfResidual(const std::vector<tREAL8> & aVResidual) const
{
   tREAL8 aSumSquare = 0;

   for (auto & aResidual : aVResidual)
       aSumSquare += Square(aResidual);

   if (mWithThr && (aSumSquare > mSig2Thrs))
      return 0.0;

   if (!mWithAtt)
      return  mWGlob;

   return  mWGlob /  (1.0 + std::pow(aSumSquare/mSig2Att,mExpS2) );
}

tREAL8  cStdWeighterResidual::SingleWOfResidual(const cPt2dr & aPt) const
{
     return SingleWOfResidual(aPt.ToStdVector());
}

std::vector<tREAL8> cStdWeighterResidual::WeightOfResidual(const tStdVect & aVResidual) const
{
   return std::vector<tREAL8>( aVResidual.size() , SingleWOfResidual(aVResidual) );
}


/* ************************************************************************ */
/*                                                                          */
/*                        cLinearWeighterResidual                           */
/*                                                                          */
/* ************************************************************************ */

cLinearWeighterResidual::cLinearWeighterResidual(tREAL8 aSGlob,tREAL8 aThr1,tREAL8 aThr2) :
    mWGlob       (1/Square(aSGlob)),
    mThrs1       (aThr1),
    mThrs2       (aThr2)
{
    if (mThrs1>mThrs2)
        std::swap(mThrs1, mThrs2);
}

tREAL8   cLinearWeighterResidual::SingleWOfResidual(tREAL8 aRes) const
{
    aRes = fabs(aRes);
    if (aRes<=mThrs1)
        return mWGlob;
    if (aRes>=mThrs2)
        return 0;
    return mWGlob*(1-(abs(aRes)-mThrs1)/(mThrs2-mThrs1));
}

std::vector<tREAL8> cLinearWeighterResidual::WeightOfResidual(const tStdVect & aVResidual) const
{
    std::vector<tREAL8> aVW;
    aVW.reserve(aVResidual.size());
    for (auto & aRes: aVResidual)
        aVW.push_back(SingleWOfResidual(aRes));
    return aVW;
}


/* ************************************************************************ */
/*                                                                          */
/*                            cMMVII_BundleAdj                              */
/*                                                                          */
/* ************************************************************************ */

cMMVII_BundleAdj::cMMVII_BundleAdj(cPhotogrammetricProject * aPhp) :
    mPhProj           (aPhp),
    mPhaseAdd         (true),
    mSys              (nullptr),
    mR8_Sys           (nullptr),
    mPatParamFrozenCalib (""),
    mPatFrozenCenter (""),
    mPatFrozenOrient (""),
    mPatFrozenClinos (""),
    //mMesGCP           (nullptr),
    //mSigmaGCP         (-1),

    mTopo             (nullptr),
    mFolderRefCam     (""),
    mSigmaTrRefCam    (-1.0),
    mSigmaRotRefCam   (-1.0),
    mPatternRef       (".*"),
    mDirRefCam        (nullptr),

    mUseGauge         (false),
    mKPoseMainGauge   (-1),
    mKPoseSecondGauge (-1),
    mKCoordSecondGauge (-1),
    mSigmaViscAngles  (-1.0),
    mSigmaViscCenter  (-1.0),
    mNbMaxIter        (-1),
    mIter             (0),
    mVerbose          (true),
    mShow_UC_UK       (false),
    mCompute_Uncert   (false),
    mRUCSUR           (nullptr),
    mVecLineAdjust    (),
    mNbCamPC          (0)
{
}


cMMVII_BundleAdj::~cMMVII_BundleAdj()
{

    ShowLVMFrozenVar();

    mSetIntervUK.SIUK_Reset();
    delete mSys;
    // delete mMesGCP;
    DeleteAllAndClear(mVTieP);
    delete mTopo;
    delete mRUCSUR;
    // DeleteAllAndClear(mGCP_UK);
    DeleteAllAndClear(mVBA_Lidar);
    DeleteAllAndClear(mVBA_LidarLidar);

    DeleteLineAdjust();
    DeleteBlockInstr();
}

void cMMVII_BundleAdj::ShowLVMFrozenVar()
{
    // surprinsinglyy cMMVII_BundleAdj can be constructed but not used in some appli,
    // this make fail the GenArgSpec
    if (mSys==nullptr)    return;

    if (mSys->NbLVMFrozen() ==0) return;

    StdOut() << " ------------  VAR W/O OBSERVATION FROZEN IN LVM -------------\n";
    size_t aIndV0 = 0;

    // Process parameters
    std::string aPatType = GetDef(mParam_UC_UK,0,std::string(".*"));  // Type selection, def=all
    std::string aPatName = GetDef(mParam_UC_UK,1,std::string(".*"));  // NameGroup selection, def=all
    std::string aPatVar =  GetDef(mParam_UC_UK,2,std::string(".*"));  // NameVar selection, def=all
    mCompute_Uncert = cStrIO<bool>::FromStr(GetDef(mParam_UC_UK,3,std::string("1")));  // Compute Uncert, def=true

    for (size_t aKObj=0 ; aKObj<  mSetIntervUK.NumberObject() ; aKObj++)
    {
        cObjWithUnkowns<tREAL8> & anObj = mSetIntervUK.KthObj(aKObj);
        cGetAdrInfoParam<tREAL8> aGIP (".*",anObj,false); // extract information
        size_t aNbVLoc = aGIP.VNames().size();

        for (size_t aKLoc=0 ; aKLoc< aNbVLoc ; aKLoc++)
        {
            size_t aKSys =  aKLoc +aIndV0;
            if (mSys->LVMFrozenUsed(aKSys))
            {
                 StdOut()  << "   * "
                           <<  " K="<< aKSys
                            <<  " Type=" << aGIP.NameType()
                             << " IdO="  << aGIP.IdObj()
                             << " Var="  <<  aGIP.VNames().at(aKLoc)
                           << "\n";
            }

        }
        aIndV0 += aNbVLoc;
    }
}

void cMMVII_BundleAdj::ShowUKNames(const std::vector<std::string> & aParam, const std::string &aSuffix, cMMVII_Appli * anAppli)
{
     // StdOut() << "=================== ShowUKNamesShowUKNames "<< aParam << " ===============\n";
     cDenseVect<tREAL8>   aVUk = mSetIntervUK.GetVUnKnowns() ;
     std::string aIdCSV = "BundleUK";
     if (!aSuffix.empty())
         aIdCSV += "_" + aSuffix;
     if (anAppli)
        anAppli->InitReportCSV(aIdCSV,"csv",false,{"Type","Group","Var","Value","Uncert"});
        // void  InitReportCSV(const std::string &anId,const std::string & aPostfix,bool IsMul,const std::vector<std::string> & aHeader={});
        // void  AddOneReportCSV(const std::string &anId,const std::vector<std::string> & VecMsg);



     for (const auto & aBBNV : mVBBNamedV)
     {
         // StdOut() << "    ************ " <<  aBBNV.mType << " : " << aBBNV.mIdObj  << "\n";
         for (size_t aKV=0 ; aKV<aBBNV.mNamesVar.size() ; aKV++)
         {
             if (aBBNV.mActivVar.at(aKV))
             {
                std::vector<std::string> aVCVS{aBBNV.mType, aBBNV.mIdObj};
                int aIndGlob = aBBNV.mIndVar0 + aKV;
                // StdOut() << "      N=" << aBBNV.mNamesVar.at(aKV)  << " V=" << aVUk(aBBNV.mIndVar0 + aKV) ;
                aVCVS.push_back(aBBNV.mNamesVar.at(aKV));
                aVCVS.push_back(ToStr(aVUk(aBBNV.mIndVar0 + aKV)));
                if (mRUCSUR)
                {
                   // StdOut()  << " UC=" <<aUC;
                   if (!mR8_Sys->VarIsFrozen(aIndGlob))
                   {
                      tREAL8 aUC = std::sqrt(mRUCSUR->UK_VarCovarEstimate(aIndGlob,aIndGlob));
                      aVCVS.push_back(ToStr(aUC));
                   }
                   else
                   {
                      aVCVS.push_back("***");
                   }
                }
                // StdOut() << "\n";
                if (anAppli)
                   anAppli->AddHeaderReportCSV(aIdCSV,aVCVS);
             }
         }
     }
     // StdOut() << "=================== ShowUKNamesShowUKNames "<< aParam << " ===============\n";
}

cPt3dr cMMVII_BundleAdj::GetGCP_UC_UK(const std::string & aGCPName) const
{
    cPt3dr aSigma = cPt3dr::Dummy();
    cDenseVect<tREAL8>   aVUk = mSetIntervUK.GetVUnKnowns() ;
    for (const auto & aBBNV : mVBBNamedV)
    {

        if (aBBNV.mIdObj == aGCPName)
        for (size_t aKV=0 ; aKV<aBBNV.mNamesVar.size() ; aKV++)
        {
                // StdOut() << "    ************ " <<  aBBNV.mType << " : " << aBBNV.mIdObj  << "\n";
            if (aBBNV.mActivVar.at(aKV))
            {
                int aIndGlob = aBBNV.mIndVar0 + aKV;
                //StdOut() << "      N=" << aBBNV.mNamesVar.at(aKV)  << " V=" << aVUk(aBBNV.mIndVar0 + aKV) << "\n";
                if (mRUCSUR)
                {
                    // StdOut()  << " UC=" <<aUC;
                    if (!mR8_Sys->VarIsFrozen(aIndGlob))
                    {
                        tREAL8 aUC = std::sqrt(mRUCSUR->UK_VarCovarEstimate(aIndGlob,aIndGlob));
                        if (aBBNV.mNamesVar.at(aKV)=="x")
                            aSigma.x() = aUC;
                        if (aBBNV.mNamesVar.at(aKV)=="y")
                            aSigma.y() = aUC;
                        if (aBBNV.mNamesVar.at(aKV)=="z")
                            aSigma.z() = aUC;
                    }
                }
            }
        }
    }
    return aSigma;
}

void cMMVII_BundleAdj::Set_UC_UK(const std::vector<std::string> & aParam)
{
     mShow_UC_UK    = true;
     mParam_UC_UK   = aParam;
}



void cMMVII_BundleAdj::AssertPhaseAdd()
{
    MMVII_INTERNAL_ASSERT_tiny(mPhaseAdd,"Mix Add and Use of unknown in cMMVII_BundleAdj");
}
void cMMVII_BundleAdj::AssertPhp()
{
    MMVII_INTERNAL_ASSERT_tiny(mPhProj,"No cPhotogrammetricProject");
}

void cMMVII_BundleAdj::AssertPhpAndPhaseAdd()
{
        AssertPhaseAdd();
        AssertPhp();
}


void cMMVII_BundleAdj::InitIteration()
{
    CompileSharedIntrinsicParams(true);
    mPhaseAdd = false;

    InitItereGCP();
    InitItereTopo();
    mR8_Sys = new cResolSysNonLinear<tREAL8>(eModeSSR::eSSR_LsqNormSparse,mSetIntervUK.GetVUnKnowns());

    mSys =  mR8_Sys;
    CompileSharedIntrinsicParams(false);

    // check if we have to add GCP UC_UK to export
    if (getGCP().getDoComputeGCP_UC_UK())
    {
        mShow_UC_UK = true;
        if (mParam_UC_UK.empty())
            mParam_UC_UK.push_back("GCP");
        else
            mParam_UC_UK.at(0)="("+mParam_UC_UK.at(0)+"|GCP)";
    }

    if (mShow_UC_UK)
    {
       size_t aIndV0 = 0;

       // Process parameters
       std::string aPatType = GetDef(mParam_UC_UK,0,std::string(".*"));  // Type selection, def=all
       std::string aPatName = GetDef(mParam_UC_UK,1,std::string(".*"));  // NameGroup selection, def=all
       std::string aPatVar =  GetDef(mParam_UC_UK,2,std::string(".*"));  // NameVar selection, def=all
       mCompute_Uncert = cStrIO<bool>::FromStr(GetDef(mParam_UC_UK,3,std::string("1")));  // Compute Uncert, def=true

       //StdOut() <<"All unknowns:\n";
       // Parse all "object"
       for (size_t aKObj=0 ; aKObj<  mSetIntervUK.NumberObject() ; aKObj++)
       {
           cObjWithUnkowns<tREAL8> & anObj = mSetIntervUK.KthObj(aKObj);

           cGetAdrInfoParam<tREAL8> aGIP (".*",anObj,false); // extract information
           cBundleBlocNamedVar aBBNV;
           aBBNV.mType = aGIP.NameType();
           aBBNV.mIdObj = aGIP.IdObj();
           aBBNV.mIndVar0 = aIndV0;
           aBBNV.mNamesVar =  aGIP.VNames();
           //StdOut() <<"    "<<aGIP.NameType()<<" "<<aGIP.IdObj()<<": ";
           //for (auto & aS:aGIP.VNames())
           //    StdOut() << aS << " ";
           //StdOut() << "\n";
           
           if (MatchRegex(aBBNV.mType,aPatType) && MatchRegex(aBBNV.mIdObj,aPatName) ) // If type and ident match
           {
               int aNbOk=0;
               for (size_t aKV=0 ; aKV<aBBNV.mNamesVar.size() ; aKV++)
               {
                   bool isOk = MatchRegex(aBBNV.mNamesVar[aKV],aPatVar);
                   aBBNV.mActivVar.push_back(isOk);
                   if (isOk)
                   {
                       aNbOk ++;
                       mIndCompUC.push_back(aKV+aIndV0);
                   }
               }
               if (aNbOk!=0)
                  mVBBNamedV.push_back(aBBNV);
           }
           aIndV0 += aBBNV.mNamesVar.size();
       }
    }
}

void cMMVII_BundleAdj::Iterate(int aNbMaxIter, tREAL8 mLVM, bool aShow_Cond)
{
    mNbMaxIter = aNbMaxIter;
    for (int aKIter=0 ; aKIter<NbMaxIter() ; aKIter++)
    {
        OneIteration(aKIter==0,mLVM,aShow_Cond);
    }
}

tREAL8   cMMVII_BundleAdj::CurLVMParam() const {return mCurLVMParam;}


void cMMVII_BundleAdj::OneIteration(bool isFirstIter, tREAL8 aLVM, bool doShowCond)
{
    mCurLVMParam = aLVM;
    MMVII_INTERNAL_ASSERT_tiny (mNbMaxIter>0,"BA: mNbMaxIter is not initialized");

    bool isLastIter =  (mIter==(mNbMaxIter-1)) ;
    // if it's first step, alloc ressources
    if (mPhaseAdd)
    {
        InitIteration();
        CheckGCPConstraints();
    }


    // ================================================
    //  [1]   Add "Hard" constraint
    // ================================================

    // if necessary, fix frozen parameters of internal calibration
    if (mPatParamFrozenCalib !="")
    {
        for (const  auto & aPtrCal : mVPCIC)
        {
            mR8_Sys->SetFrozenFromPat(*aPtrCal,mPatParamFrozenCalib,true);
        }
    }
    for (const auto & aPat : mPatternFreeCalib)
    {
        for (const  auto & aPtrCal : mVPCIC)
        {
            if (MatchRegex(aPtrCal->Name(),aPat.PatCal))
            {
               mR8_Sys->SetFrozenFromPat(*aPtrCal,aPat.PatParam,false,aPat.Weight);
            }
        }
    }

    if (mUseGauge)
    {
        cSensorCamPC * aMainCamGauge   = mVSCPC.at(mKPoseMainGauge);
        cSensorCamPC * aSecondCamGauge = mVSCPC.at(mKPoseSecondGauge);

        mR8_Sys->SetFrozenVarCurVal(*aMainCamGauge,aMainCamGauge->Center());
        mR8_Sys->SetFrozenVarCurVal(*aMainCamGauge,aMainCamGauge->Omega());

        tREAL8 & aCoord = aSecondCamGauge->Center()[mKCoordSecondGauge];
        mR8_Sys->SetFrozenVarCurVal(*aSecondCamGauge,aCoord);

    }

    // if necessary, fix frozen centers of external calibration
    if (mPatFrozenCenter !="")
    {
        tNameSelector aSel =   AllocRegex(mPatFrozenCenter);
        int nbMatches = 0;
        for (const auto & aPtrCam : mVSCPC)
        {
            if ((aPtrCam != nullptr)  && aSel.Match(aPtrCam->NameImage()))
            {
                mR8_Sys->SetFrozenVarCurVal(*aPtrCam,aPtrCam->Center());
                nbMatches++;
            }
        }
        if (mVerbose && isFirstIter)
            StdOut() << "Frozen centers: " << nbMatches << ".\n";
    }
   
    // if necessary, fix frozen orientation of external calibration
    if (mPatFrozenOrient !="")
    {
        tNameSelector aSel =   AllocRegex(mPatFrozenOrient);
        int nbMatches = 0;
        for (const auto & aPtrCam : mVSCPC)
        {
            if ((aPtrCam != nullptr)  && aSel.Match(aPtrCam->NameImage()))
            {
                mR8_Sys->SetFrozenVarCurVal(*aPtrCam,aPtrCam->Omega());
                nbMatches++;
            }
        }
        if (mVerbose && isFirstIter)
            StdOut() << "Frozen orients: " << nbMatches << ".\n";
    }
    // if necessary fix hard cosntraint onf Gauge of Rigid-Block of instrument
    SetHardGaugeBlockInstr();


    if (mTopo) // TOPO
    {
        mTopo->SetFrozenAndSharedVars(*mR8_Sys);
    }

    // ================================================
    //  [2]   Add "Soft" constraint
    // ================================================

    // if necessary, add some "viscosity" on poses
    AddPoseViscosity();

    // Add constriant betweenn reference and pose
    AddConstrainteRefPose();


    // ================================================
    //  [3]   Add compensation measures
    // ================================================


    OneItere_GCP();   // add GCP informations
    OneItere_TieP();  // ad tie-points information


                      //


    for (const auto & aLidarPh : mVBA_Lidar )
        aLidarPh->AddObs();

    for (const auto & aLidarLidar : mVBA_LidarLidar )
        aLidarLidar->AddObs();

    // Add observation for line adjustment
    IterAdjustOnLine();
    //  Add observation for block of instrument
    IterOneBlockInstr();

    // add topo last to compute sigma0 will all constraints
    if (mTopo) // TOPO
    {
        mTopo->AddTopoEquations(*mR8_Sys);
#ifdef VERBOSE_TOPO
        mTopo->print();
#endif
        if (mVerbose)
            mTopo->printObs(false);
    }

    if (mCompute_Uncert && isLastIter)
    {
// StdOut() <<  "mCompute_UncertmCompute_UncertmCompute_Uncert--------------------------------\n";
        mRUCSUR = new cResult_UC_SUR<tREAL8>(false,false,mIndCompUC);
    }

    mR8_Sys->SysLinear(); // to force SetPhaseEq(); to init constraints if no AddObs

    const auto & aVectSol = mR8_Sys->SolveUpdateReset(aLVM,{},{mRUCSUR},doShowCond);
    mSetIntervUK.SetVUnKnowns(aVectSol);

    mIter++;
    if(mVerbose)
    {
        StdOut() << "---------------------- "
                 << " End Iter" << mIter
                  << " StdDevLast=" << std::sqrt(mR8_Sys->VarLastSol())
                  << " StdDevCur=" << std::sqrt(mR8_Sys->VarCurSol())
                  //<< " VarLast=" << mR8_Sys->VarLastSol()
                 //<< " ?VarCur?=" << mR8_Sys->VarCurSol()
                 << " ---------------" << std::endl;
    }

}


void cMMVII_BundleAdj::AddCalib(cPerspCamIntrCalib * aCalib)
{
    AssertPhaseAdd();
    if (aCalib==nullptr)  return;
    if (! aCalib->UkIsInit())
    {
          mVPCIC.push_back(aCalib);
          mSetIntervUK.AddOneObj(aCalib);
    }
}

void  cMMVII_BundleAdj::AddBenchSensor(cSensorCamPC * aSI){
    mVSCPC.push_back(aSI);
    AddSensor(aSI);
}

void cMMVII_BundleAdj::AddSensor(cSensorImage* aSI)
{
    AssertPhaseAdd();
    MMVII_INTERNAL_ASSERT_tiny (!aSI->UkIsInit(),"Multiple add of cam : " + aSI->NameImage());
    mSetIntervUK.AddOneObj(aSI);
    mVSIm.push_back(aSI);

    aSI->SetAndGetEqColinearity(true,10,true);  // WithDer, SzBuf, ReUse
}


void cMMVII_BundleAdj::AddCamPC(cSensorCamPC * aCamPC)
{
    mNbCamPC++;
    if (aCamPC->DoAddCalibToUk())
        AddCalib(aCamPC->InternalCalib());
}

void cMMVII_BundleAdj::AddReferencePoses(const std::string & aOri, tREAL8 aSigmaTr, tREAL8 aSigmaRot, const std::string & aPatApply)
{
     MMVII_INTERNAL_ASSERT_tiny(mVSCPC.empty(),"Must Add Ref Pose before any cam");
     AssertPhpAndPhaseAdd();

     mFolderRefCam = aOri;
     mDirRefCam  = mPhProj->NewDPIn(eTA2007::Orient,mFolderRefCam);

     mSigmaTrRefCam = aSigmaTr;
     mSigmaRotRefCam = aSigmaRot;
     mPatternRef = aPatApply;
}


void  cMMVII_BundleAdj::AddCam(const std::string & aNameIm)
{
    AssertPhpAndPhaseAdd();

    cSensorImage * aNewS;
    cSensorCamPC * aSPC;

    mPhProj->ReadSensor(aNameIm,aNewS,aSPC,true,false);  // false -> NoSVP
    AddSensor(aNewS);

    mVSCPC.push_back(aSPC);  // eventually nullptr, for example with push-broom
    if (aSPC)
       AddCamPC(aSPC);

    //  process the reference cameras
    if (mDirRefCam)
    {
        // if PC Cam dont exist push 0 to be coherent between "mVCamRefPose" and   "mVSCPC"
        if (aSPC==nullptr)
        {
            mVCamRefPoses.push_back(aSPC);
        }
        else
        {
            const std::string & aNameImage = aSPC->NameImage();
            cSensorCamPC * aCamRef  = mPhProj->ReadCamPC(*mDirRefCam,aNameImage,true);
            mVCamRefPoses.push_back(aCamRef);
        }
    }
}
const std::vector<cSensorImage *> &  cMMVII_BundleAdj::VSIm() const  {return mVSIm;}
const std::vector<cSensorCamPC *> &  cMMVII_BundleAdj::VSCPC() const {return mVSCPC;}
cResolSysNonLinear<tREAL8> *  cMMVII_BundleAdj::Sys() {return mR8_Sys;}
cPhotogrammetricProject  & cMMVII_BundleAdj::PhProj() {return *mPhProj;}


cSetInterUK_MultipeObj<tREAL8> &   cMMVII_BundleAdj::SetIntervUK() {return mSetIntervUK;}
int cMMVII_BundleAdj::NbCamPC() const {return mNbCamPC;}

    /* ---------------------------------------- */
    /*            Frozen/Shared                 */
    /* ---------------------------------------- */

void cMMVII_BundleAdj::SetParamFrozenCalib(const std::string & aPattern)
{
    mPatParamFrozenCalib = aPattern;
}

void cMMVII_BundleAdj::SetParamFreeCalib(const std::vector<cFreeCalibPattern> & aVPat)
{
    mPatternFreeCalib = aVPat;
}

void cMMVII_BundleAdj::SetFrozenCenters(const std::string & aPattern)
{
    mPatFrozenCenter = aPattern;
}

void cMMVII_BundleAdj::SetFrozenOrients(const std::string & aPattern)
{
    mPatFrozenOrient = aPattern;
}

void cMMVII_BundleAdj::SetFrozenClinos(const std::string & aPattern)
{
    mPatFrozenClinos = aPattern;
}

void cMMVII_BundleAdj::SetSharedIntrinsicParams(const std::vector<cSharedIPParam> & aVParams)
{
    mVPatShared = aVParams;
}

typedef std::tuple<int,std::string,tREAL8 *> tISRP;

void cMMVII_BundleAdj::CompileSharedIntrinsicParams(bool ForAvg)
{
    bool  Show = ForAvg;

    // Parse the pair Pattern Name Cam / Pattern Name Params
    for (const auto & aPatShared : mVPatShared)
    {
        std::map<std::string,std::vector<int>> aMapSharedIndexes; // store the shared index of a given param name
        std::map<std::string,std::vector<std::string>> aMapNames; // store the sharing as name, for show
        std::map<std::string,std::vector<tISRP>> aMapValues; // store the sharing of adress for averaging
        // Parse the calib and select those which name match the pattern name cam
        for (auto  aPtrCal : mVPCIC)
        {
            if (MatchRegex(aPtrCal->Name(),aPatShared.PatCal))
            {
                // Extract information on parameter macthing the pattern of params
                cGetAdrInfoParam<tREAL8>  aGIP(aPatShared.PatParam,*aPtrCal,false);
                for (size_t aKParam=0 ; aKParam<aGIP.VAdrs().size() ; aKParam++)
                {
                    tREAL8 * aAdr                 = aGIP.VAdrs().at(aKParam);
                    const std::string & aNameP    = aGIP.VNames().at(aKParam);
                    cObjWithUnkowns<tREAL8>* aObj = aGIP.VObjs().at(aKParam);

                    size_t  aNum = aObj->IndOfVal(aAdr);
                    aMapSharedIndexes[aNameP].push_back(aNum);
                    aMapNames[aNameP].push_back(aPtrCal->Name());

                    aMapValues[aNameP].push_back(tISRP(aNum,aPtrCal->Name(),aAdr));
                }
            }
        }
        if (Show)
        {
                StdOut()  << "=========== Shared params for"
                          << " PatCal={" << aPatShared.PatCal << "}"
                          << " PatParam={" << aPatShared.PatParam << "}"
                          << " ============ " << std::endl;
        }
        for (const auto & [aNamePar,aVTuple] : aMapValues)
        {
            std::vector<int>  aVIndEqui;
            tREAL8  aSum = 0.0;
            for (const auto & [aNum,aNameCam,anAdr] : aVTuple)
            {
                aSum += *anAdr;
                aVIndEqui.push_back(aNum);
            }
            aSum /= aVTuple.size();
            if (Show)
               StdOut()  <<  "   * "  << aNamePar  << " : " << aSum << std::endl;
            if (ForAvg)
            {
                for (const auto & [aNum,aNameCam,anAdr] : aVTuple)
                {
                    //*anAdr = aSum;  => dont undertand why it apparently slow down the convergence ??
                    StdOut()  <<  "      - "  << aNameCam <<  " : " << *anAdr << std::endl;
                }
            }
            else
            {
               mSys->SetShared(aVIndEqui);
            }
        }
/*
        if (ForAvg)
        {
           if (Show)
           {
                StdOut()  << "=========== Shared params for"
                          << " PatCal={" << aPatShared.PatCal << "}"
                          << " PatParam={" << aPatShared.PatParam << "}"
                          << " ============ " << std::endl;
                for (const auto & [aNamePar,aVNameCam] : aMapNames)
                {
                      StdOut()  <<  "   * "  << aNamePar << std::endl;
                      for (const auto &  aNameCam : aVNameCam )
                          StdOut()  <<  "      - "  << aNameCam << std::endl;
                }
                StdOut()  << "==========================================================" << std::endl;
           }
        }
        else
        {
             for (const auto & [aNamePar,aVIndexes] : aMapSharedIndexes)
             {
                 mSys->SetShared(aVIndexes);
             }
        }
*/
    }
}

int cMMVII_BundleAdj::IndexOfPCPose(const std::string &aNameIm,bool SVP ) const
{
    for (size_t aKP=0 ; aKP<mVSCPC.size() ; aKP++)
    {
        if (mVSCPC.at(aKP) && (mVSCPC.at(aKP)->NameImage()==aNameIm))
            return aKP;
    }
    MMVII_INTERNAL_ASSERT_always(SVP,"IndexOfPCPose, cannot get name :"+aNameIm);
    return -1;
}

void cMMVII_BundleAdj::SetGaugeRelPause(const cParamFixGauge & aParam)
{

    cWhichMax<cPt3di,tREAL8> aWMaxInd;

    //  an empty field means "not fixed" : the system then chooses it
    const std::string & aN1 = aParam.MainIm;
    const std::string & aN2 = aParam.SecIm;

    bool aN1Fix =    ! aN1.empty();
    bool aN2Fix =    ! aN2.empty();
    bool aCoordFix = ! aParam.Coord.empty();

    std::vector<std::string> aVCoord{"x","y","z"};
    int aKCoord = -1;
    if (aCoordFix)
    {
        auto anIter = std::find(aVCoord.begin(),aVCoord.end(),aParam.Coord);
        MMVII_INTERNAL_ASSERT_always(anIter!=aVCoord.end(),"SetGaugeRelPause bad coord");
        aKCoord = anIter - aVCoord.begin();
    }

    int aKCoordBegin = aCoordFix      ?  aKCoord    : 0 ;
    int aKCoordEnd   = aCoordFix      ? (aKCoord+1) : 3 ;



    for (size_t aKP1=0 ; aKP1<mVSCPC.size() ; aKP1++)
    {
        cSensorCamPC * aCam1 = mVSCPC.at(aKP1);
        if (aCam1 && ((!aN1Fix) || (aCam1->NameImage()==aN1))  )
        {
           cPt3dr aC1 = aCam1->Center();
           int aBeginKP2 = aN1Fix ? 0 : (aKP1+1);
           for (size_t aKP2=aBeginKP2 ; aKP2<mVSCPC.size() ; aKP2++)
           {
                cSensorCamPC * aCam2 = mVSCPC.at(aKP2);
                if (aCam2 && ((!aN2Fix) || (aCam2->NameImage()==aN2)) )
                {
                    cPt3dr aV12 = aCam2->Center() - aC1;
                    for (int aKCoord=aKCoordBegin; aKCoord<aKCoordEnd ; aKCoord++)
                    {
                        aWMaxInd.Add(cPt3di(aKP1,aKP2,aKCoord),std::abs(aV12[aKCoord]));
                    }
                }
           }
        }
    }

    MMVII_INTERNAL_ASSERT_always(aWMaxInd.IsInit(),"Cannot compute relative gauge");

    cPt3di aIndMax = aWMaxInd.IndexExtre();

    if (1)
    {
        StdOut() << "GAUGE SET "
                 << " Im1=" << mVSCPC.at(aIndMax.x())->NameImage()
                 << " Im2=" << mVSCPC.at(aIndMax.y())->NameImage()
                 << " Coord=" << aVCoord.at(aIndMax.z())
                 << " Dist=" << aWMaxInd.ValExtre()
                 << "\n";


    }
    mUseGauge = true;
    mKPoseMainGauge = aIndMax.x();
    mKPoseSecondGauge = aIndMax.y();
    mKCoordSecondGauge = aIndMax.z();


//    aMaxInd
}


//void SetGaugeRelPause(int aKPoseMain,int aKposeSec,int aKCoord);


bool cMMVII_BundleAdj::CheckGCPConstraints() const
{
    std::string aNames;

    for (const auto & aMesGCP : mGCP.getMesGCP().MesGCP())
    {
        if (aMesGCP.isFree())
        {

            int aNbImObs =  mGCP.getMesGCP().GetNbImMesForPoint(aMesGCP.mNamePt);
            int aNbTopoElementObs = 0;
            if (mTopo)
            {
                for (auto & obs: mTopo->GetObsPoint(aMesGCP.mNamePt))
                {
                    aNbTopoElementObs += obs->getMeasures().size();
                }
            }
            if (aNbImObs+aNbTopoElementObs<3)
                aNames += aMesGCP.mNamePt + " ";
        }
    }
    if (aNames.size())
        MMVII_USER_WARNING("Not enough observations for points: "+aNames);

    return true;
}

    /* ---------------------------------------- */
    /*            AddViscosity                  */
    /* ---------------------------------------- */

void cMMVII_BundleAdj::AddPoseViscosity()
{
     //  parse all centra
     for (auto aPcPtr : mVSCPC)
     {
         if (aPcPtr!=nullptr)
         {
            if (mSigmaViscCenter>0)
            {
               mR8_Sys->AddEqFixCurVar(*aPcPtr,aPcPtr->Center(),Square(1/mSigmaViscCenter));
            }
            if (mSigmaViscAngles>0)
            {
               mR8_Sys->AddEqFixCurVar(*aPcPtr,aPcPtr->Omega(),Square(1/mSigmaViscAngles));
            }
         }
     }
}


void cMMVII_BundleAdj::SetViscosity(const tREAL8& aViscTr,const tREAL8& aViscAngle)
{
    mSigmaViscCenter = aViscTr;
    mSigmaViscAngles = aViscAngle;
}

    /* ---------------------------------------- */
    /*             Reference Pose               */
    /* ---------------------------------------- */

void cMMVII_BundleAdj::AddConstrainteRefPose()
{
   if (!mDirRefCam)
      return;

   for (size_t aKC=0 ; aKC<mVSCPC.size() ; aKC++)
   {
        cSensorCamPC * aCam = mVSCPC[aKC];
        cSensorCamPC * aCamRef =  mVCamRefPoses[aKC];
        if ((aCam!=nullptr) && (aCamRef!=nullptr))
           AddConstrainteRefPose(*aCam,*aCamRef);
   }
}

void cMMVII_BundleAdj::AddConstrainteRefPose(cSensorCamPC & aCam,cSensorCamPC & aCamRef)
{
     if (! MatchRegex(aCam.NameImage(),mPatternRef))
        return;
     // mR8_Sys
     if (mSigmaTrRefCam > 0)
        mR8_Sys->AddEqFixNewVal(aCam,aCam.Center(),aCamRef.Center(),Square(1/mSigmaTrRefCam));
     
     if (mSigmaRotRefCam>0)
     {
         cPt3dr aWTarget = aCam.Pose_WU().ValAxiatorFixRot(aCamRef.Pose().Rot());
         mR8_Sys->AddEqFixNewVal(aCam,aCam.Omega(),aWTarget,Square(1/mSigmaRotRefCam));
     }

}

    /* ---------------------------------------- */
    /*            Rigid Bloc                    */
    /* ---------------------------------------- */




/* ---------------------------------------- */
/*                 Lidar                    */
/* ---------------------------------------- */

void cMMVII_BundleAdj::Add1AdjLidarPhotogra(eImatchCrit aMode, const std::string & aPlyFile, double aSigma,
                                            const std::vector<std::string> & aInterp, bool aPertubate, int aNbPtsPerPatch)
{
    mVBA_Lidar.push_back(new cBA_LidarPhotograTri(mPhProj, *this, aMode, aPlyFile, aSigma, aInterp, aPertubate, aNbPtsPerPatch));
}

void cMMVII_BundleAdj::Add1AdjLidarPhoto(eImatchCrit aMode, const std::string & aPatScan, double aSigma,
                                         const std::vector<std::string> & aInterp, double aScaleInit, double aScaleFinal,
                                         double aThreshold, int aNbPtsPerPatch)
{
    mVBA_Lidar.push_back(new cBA_LidarPhotograRaster(mPhProj, *this, aMode, aPatScan, aSigma, aInterp,
                                                     aScaleInit, aScaleFinal, aThreshold, aNbPtsPerPatch));
}

bool cMMVII_BundleAdj::AddStaticLidar(cStaticLidar* aStaticLidar)
{
    AssertPhpAndPhaseAdd();

    if (mMapTSL.count(aStaticLidar->NameImage())==0)
    {
        aStaticLidar->ReadRasters(mPhProj->DirStaticLidarRasters());
        MMVII_INTERNAL_ASSERT_tiny (mMapTSL.count(aStaticLidar->NameImage())==0,"Multiple add of TSL : " + aStaticLidar->NameImage());
        mMapTSL[aStaticLidar->NameImage()] = aStaticLidar;
    }
    return true;
}

const std::unordered_map<std::string, cStaticLidar *> &cMMVII_BundleAdj::MapTSL() const {return mMapTSL;}

void cMMVII_BundleAdj::Add1AdjLidarLidar(const std::string & aPatScan, double aSigma, double aThresholdInit,
                                         double aThresholdFinal, double aNormalTolDeg, const std::vector<std::string> & aInterp)
{
    mVBA_LidarLidar.push_back(new cBA_LidarLidarRaster(mPhProj, *this, aPatScan, aSigma, aThresholdInit,
                                                       aThresholdFinal, aNormalTolDeg, aInterp));
}


/* ---------------------------------------- */
/*                 Topo                     */
/* ---------------------------------------- */

void cMMVII_BundleAdj::SaveTopo()
{
    if (mTopo && mPhProj->DPTopoMes().DirOutIsInit())
    {
        mPhProj->SaveTopoMes(*mTopo);
    }
}


void cMMVII_BundleAdj::AddTopo() // TOPO
{
    mTopo = new cBA_Topo(mPhProj, &mGCP);
    mTopo->AddPointsFromDataToGCP( mPhProj);
}

}; // MMVII

