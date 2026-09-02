#include "BundleAdjustment.h"  // also pulls MMVII_Tpl_ElemStrToVal.h (for ARG2007_STRUCT_FIELDS)

/**
   \file cAppliBundAdj.cpp

*/
#include "MMVII_Topo.h"

/*
    Track info on bundle adj/push broom in V1 :

  * mm3d Campari

       * [Name=PdsGBRot] REAL :: {Weighting of the global rotation constraint (Generic bundle Def=0.002)}
       * [Name=PdsGBId] REAL :: {Weighting of the global deformation constraint (Generic bundle Def=0.0)}
       * [Name=PdsGBIter] REAL :: {Weighting of the change of the global rotation constraint between iterations (Generic bundle Def=1e-6)}

  *  micmac/src/uti_phgrm/CPP_Campari.cpp [POS=0486,0016]
    
          Apero ... "Apero-Compense.xml" ...
             +  std::string(" +PdsGBRot=") + ToString(aPdsGBRot) + " "

  * micmac/include/XML_MicMac/Apero-Compense.xml

       <ContrCamGenInc>
             <PatternApply> .*  </PatternApply>
             <PdsAttachToId>   ${PdsGBId}     </PdsAttachToId>
             <PdsAttachToLast> ${PdsGBIter}    </PdsAttachToLast>
             <PdsAttachRGLob>  ${PdsGBRot}    </PdsAttachRGLob>
        </ContrCamGenInc>


  *  micmac/src/uti_phgrm/Apero/cPosePolynGenCam.cpp

       cPolynBGC3M2D_Formelle * aPF = aGPC->PolyF();

       if (aCCIG.PdsAttachRGLob().IsInit())
          aPF->AddEqRotGlob(aCCIG.PdsAttachRGLob().Val()*aGPC->SomPM());


   * micmac/src/uti_phgrm/Apero/cGenPoseCam.cpp
   * micmac/src/uti_phgrm/Apero/BundleGen.h

*/


namespace MMVII
{

struct cGCP3D
{
    std::string GCPDir;
    double Sigma;
    std::string OutDir="";
    double ExportSigma = 0;

    ARG2007_STRUCT_FIELDS (
        GCPDir,FieldSem({{eTA2007::AddCom,"World coord GCP dir"},eTA2007::ObjCoordWorld}),
        Sigma,FieldSem({{eTA2007::AddCom,"Sigma factor SG=0 fix, SG<0 schurr elim, SG>0"},eTA2007::ObjCoordWorld}),
        OutDir,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Optional output dir"}}),
        ExportSigma,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Optional compensated sigma"}})
    )
};

struct cGCP2D
{
    std::string GCPDir;
    double Sigma;
    double SigmaAttenuation=-1;
    double Threshold=-1;
    double Exponent=1;

    ARG2007_STRUCT_FIELDS (
        GCPDir,FieldSem({{eTA2007::AddCom,"Image coords GCP dir"},eTA2007::ObjMesInstr}),
        Sigma,FieldSem({eTA2007::AddCom,"Sigma factor"}),
        SigmaAttenuation,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Optional Sigma attenuation"}}),
        Threshold,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Optional Sigma threshold"}}),
        Exponent,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Optional Sigma exponent"}})
        )
};

struct cLidarPhotograParam
{
    eImatchCrit Mode;
    std::string PlyFile;
    double Sigma;
    std::vector<std::string> Interp = {"Tabul","1000","SinCApod","10","10"};
    bool Pertubate = false;
    int NbPtsPerPatch = 49;
    ARG2007_STRUCT_FIELDS (
        Mode,FieldSem({eTA2007::AddCom,"Similarity criterion used to compare Lidar and image"}),
        PlyFile,FieldSem({eTA2007::FileCloud,{eTA2007::AddCom,"Lidar point cloud"}}),
        Sigma,FieldSem({eTA2007::AddCom,"Sigma factor"}),
        Interp,FieldSem({eTA2007::HDV,eTA2007::Interpol,{eTA2007::AddCom,"Interpolator used to sample images"}}),
        Pertubate,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Perturbate radiometry (simulation/test)"}}),
        NbPtsPerPatch,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Approximate number of points per patch"}})
        )
};

struct cLidarPhotoParam
{
    eImatchCrit Mode;
    std::string PatScan;
    double Sigma;
    std::vector<std::string> Interp = {"Tabul","1000","SinCApod","10","10"};
    double ScaleInit = 1;
    double ScaleFinal = 1;
    double Threshold = -1;
    int NbPtsPerPatch = 49;
    ARG2007_STRUCT_FIELDS (
        Mode,FieldSem({eTA2007::AddCom,"Similarity criterion used to compare Lidar and image"}),
        PatScan,FieldSem({eTA2007::AddCom,"Pattern of scan names to use"}),
        Sigma,FieldSem({eTA2007::AddCom,"Sigma factor"}),
        Interp,FieldSem({eTA2007::HDV,eTA2007::Interpol,{eTA2007::AddCom,"Interpolator used to sample images"}}),
        ScaleInit,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Interpolator scale at first iteration"}}),
        ScaleFinal,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Interpolator scale at last iteration"}}),
        Threshold,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Distance threshold for hidden points, <0 = no threshold"}}),
        NbPtsPerPatch,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Approximate number of points per patch"}})
        )
};

struct cLidarLidarParam
{
    std::string PatScan;
    double Sigma;
    double ThresholdInit = 1;
    double ThresholdFinal = 0.1;
    double NormalTolDeg = 15;
    std::vector<std::string> Interp = {"Linear"};
    ARG2007_STRUCT_FIELDS (
        PatScan,FieldSem({eTA2007::AddCom,"Pattern of scan names to use"}),
        Sigma,FieldSem({eTA2007::AddCom,"Sigma factor"}),
        ThresholdInit,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Distance threshold at first iteration, <0 = infinite"}}),
        ThresholdFinal,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Distance threshold at last iteration, <0 = infinite"}}),
        NormalTolDeg,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Max normal angle tolerance (degrees)"}}),
        Interp,FieldSem({eTA2007::HDV,eTA2007::Interpol,{eTA2007::AddCom,"Interpolator used to sample scans"}})
        )
};

struct cAddTieP
{
    std::string Folder;
    double Sigma;
    double SigmaAttenuation=-1;
    double Threshold=-1;
    double Exponent=1;
    ARG2007_STRUCT_FIELDS (
        Folder,FieldSem({{eTA2007::AddCom,"Tie points folder"},eTA2007::MulTieP}),
        Sigma,FieldSem({eTA2007::AddCom,"Sigma factor"}),
        SigmaAttenuation,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Optional Sigma attenuation"}}),
        Threshold,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Optional Sigma threshold"}}),
        Exponent,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Optional Sigma exponent"}})
        )
};

struct cTiePWeightParam
{
    double Sigma;
    double SigmaAttenuation=-1;
    double Threshold=-1;
    double Exponent=1;
    ARG2007_STRUCT_FIELDS (
        Sigma,FieldSem({eTA2007::AddCom,"Sigma factor"}),
        SigmaAttenuation,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Optional Sigma attenuation"}}),
        Threshold,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Optional Sigma threshold"}}),
        Exponent,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Optional Sigma exponent"}})
        )
};

struct cFreeCalibParam
{
    std::string PatCal;
    std::string PatParam;
    double Weight=-1;
    ARG2007_STRUCT_FIELDS (
        PatCal,FieldSem({eTA2007::AddCom,"Pattern on calibration name"}),
        PatParam,FieldSem({eTA2007::AddCom,"Pattern on parameter name"}),
        Weight,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Optional weight, <0 = hard constraint"}})
        )
};

struct cRefOriParam
{
    std::string Ori;
    double SigmaTr;
    double SigmaRot=-1;
    std::string PatApply=".*";
    ARG2007_STRUCT_FIELDS (
        Ori,FieldSem({{eTA2007::AddCom,"Reference orientation folder"},eTA2007::Orient}),
        SigmaTr,FieldSem({eTA2007::AddCom,"Sigma on translation"}),
        SigmaRot,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Optional sigma on rotation"}}),
        PatApply,FieldSem({eTA2007::HDV,{eTA2007::AddCom,"Pattern of images to apply constraint"}})
        )
};

struct cAdjLine3DParam
{
    std::string Folder;
    double SigmaIm;
    int NbPtsSampl;
    ARG2007_STRUCT_FIELDS (
        Folder,FieldSem({eTA2007::AddCom,"Folder for 3D lines"}),
        SigmaIm,FieldSem({eTA2007::AddCom,"Sigma in image"}),
        NbPtsSampl,FieldSem({eTA2007::AddCom,"Number of sampled points"})
        )
};

struct cPoseViscParam
{
    double SigmaCenter;
    double SigmaRot;
    ARG2007_STRUCT_FIELDS (
        SigmaCenter,FieldSem({eTA2007::AddCom,"Sigma viscosity on center"}),
        SigmaRot,FieldSem({eTA2007::AddCom,"Sigma viscosity on rotation"})
        )
};

   /* ********************************************************** */
   /*                                                            */
   /*                 cAppliBundlAdj                             */
   /*                                                            */
   /* ********************************************************** */

class cAppliBundlAdj : public cMMVII_Appli
{
     public :
        cAppliBundlAdj(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec);
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;
     private :
        bool AcceptEmptySet(int aK) const override {return ((aK==0)&&(mSpecImIn==MMVII_NONE));}

        void  AddOneSetGCP3D(const std::string & aFolderIn, const std::string &aFolderOut, tREAL8 aWFactor, bool aDoExportSigmas); // aFolderOut="" if no out
        void  AddOneSetGCP2D(const cGCP2D & aGCP2D);
        void  AddOneSetTieP(const std::string & aFolder, double aSigma, double aSigmaAttenuation, double aThreshold, double aExponent);

        std::string               mSpecImIn;

       // std::string               mDataDir;  /// Default Data dir for all

        cPhotogrammetricProject   mPhProj;
        cMMVII_BundleAdj          mBA;

        std::vector<cGCP3D>       mGCP3D; // gcp ground coords with sigma factor and optional output dir
        std::vector<cGCP2D>       mGCP2D; // gcp image coords with weight
        std::string               mGCPFilter;  // pattern to filter names of GCP
        std::string               mGCPFilterAdd;  // pattern to filter GCP by additional info
        cTiePWeightParam          mTiePWeight;
        std::vector<cAddTieP>     mAddTieP; // In case there is multiple GCP Set
        std::vector<int>          mTiepShowPerMil;
        cRefOriParam              mParamRefOri;  // Force Poses to be +- equals to this reference

        std::vector<cLidarPhotograParam>  mParamLidarPhgr; // parameters for photogra/lidar adj via triangulation
        std::vector<cLidarPhotoParam>     mParamLidarPhoto; // parameters for photogra/lidar adj via rasterization
        std::vector<cLidarLidarParam>     mParamLidarLidar; // parameters for lidar-lidar adj

        int                       mNbIter;

        std::string               mPatParamFrozCalib;
        std::vector<cFreeCalibParam>  mVVParFreeCalib;
        cParamFixGauge            mParamGaugeRel;
        std::string               mPatFrosenCenters;
        std::string               mPatFrosenOrient;
        std::string               mPatFrosenClino;
        cPoseViscParam            mViscPose;
        tREAL8                    mLVM;  ///< Levenberk Markard
        bool                      mMeasureAdded ;
        std::vector<cSharedIPParam>  mVSharedIP;  ///< Vector for shared intrinsic param
        bool                      mShow_Cond; ///< compute and show system condition number
        std::vector<std::string>  mParamShow_UK_UC;
        std::string               mPostFixReport;
        cAdjLine3DParam           mParamLine;
        cParamBOI                 mParamBOI;  //< Param for bloc of instrum
        cParamClinoBOI            mParamBOIClino;  //< Param for clino of bloc of instrum
         /// For tuning other options, like LVM, we may want to omit the check on measure
         bool                     mCheckMeasureAdded;
};
cAppliBundlAdj::cAppliBundlAdj(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
   cMMVII_Appli    (aVArgs,aSpec),
  // mDataDir        ("Std"),
   mPhProj            (*this),
   mBA                (&mPhProj),
   mGCPFilter         (""),
   mGCPFilterAdd      (""),
   mTiepShowPerMil    {500,750},
   mNbIter            (10),
   mLVM               (0.0),
   mMeasureAdded      (false),
   mShow_Cond         (false),
   mCheckMeasureAdded (true)
{
}

cCollecSpecArg2007 & cAppliBundlAdj::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return anArgObl
              << Arg2007(mSpecImIn,"Pattern/file for images",{{eTA2007::MPatFile,"0"},{eTA2007::FileDirProj}})
              <<  mPhProj.DPOrient().ArgDirInMand()
              <<  mPhProj.DPOrient().ArgDirOutMand()
           ;
}


cCollecSpecArg2007 & cAppliBundlAdj::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    return
          anArgOpt
     // << AOpt2007(mDataDir,"DataDir","Default data directories ",{eTA2007::HDV})
      << AOpt2007(mParamRefOri,"RefOri","Reference orientation")
      << AOpt2007(mPostFixReport,NameParamPostFixReport(),CommentParamPostFixReport())
      << AOpt2007(mParamLine,"AdjLine3D","Parameter for line Adjustment")

      << cHeaderSectionArg("Topo")
      << mPhProj.DPTopoMes().ArgDirInOpt("TopoDirIn","Dir for Topo measures") //  TOPO
      << mPhProj.DPTopoMes().ArgDirOutOpt("TopoDirOut","Dir for Topo measures output") //  TOPO

      << cHeaderSectionArg("Clino")
      << mPhProj.DPMeasuresClino().ArgDirInOpt()

      << cHeaderSectionArg("GCPs")
      << AOpt2007 ( mGCP3D, "GCP3D", "GCP ground",{eTA2007::CanRepeat})
      << AOpt2007 ( mGCP2D, "GCP2D", "GCP image",{eTA2007::CanRepeat})
      << AOpt2007(mGCPFilter,"GCPFilter","Pattern to filter GCP by name")
      << AOpt2007(mGCPFilterAdd,"GCPFilterAdd","Pattern to filter GCP by additional info")

      << cHeaderSectionArg("Tie points")
      << mPhProj.DPMulTieP().ArgDirInOpt("TPDir","Dir for Tie Points if != DataDir")
      << AOpt2007(mTiePWeight,"TiePWeight","Tie point weighting")
      << AOpt2007(mAddTieP,"AddTieP","For additional TieP",{eTA2007::CanRepeat})
      << AOpt2007(mTiepShowPerMil,"TiePShowPerMil","Per/1000 for printing residual, for ex [500] -> show median ",{eTA2007::HDV})


      << cHeaderSectionArg("Lidar")
      << AOpt2007 ( mParamLidarPhgr,  "LidarPhotogra", "Lidar/Phgr adj via triangulation",  {eTA2007::CanRepeat})
      << AOpt2007 ( mParamLidarPhoto, "LidarPhoto",    "Lidar/Phgr adj via rasterisation",  {eTA2007::CanRepeat})
      << AOpt2007 ( mParamLidarLidar, "LidarLidar",    "Lidar/Lidar adj via rasterisation", {eTA2007::CanRepeat})

      << cHeaderSectionArg("Freeze")
      << AOpt2007(mPatParamFrozCalib,"PPFzCal","Pattern for freezing internal calibration parameters")
      << AOpt2007(mVVParFreeCalib,"PPFreeCal","Pattern for free internal calibration parameters",{eTA2007::CanRepeat})
      << AOpt2007(mPatFrosenCenters,"PatFzCenters","Pattern of images for freezing center of poses")
      << AOpt2007(mPatFrosenOrient,"PatFzOrient","Pattern of images for freezing orientation of poses")
      << AOpt2007(mPatFrosenClino,"PatFzClino","Pattern of clinometers for freezing boresight")
      << AOpt2007(mParamGaugeRel,"FixGauge","Param for gauge in pure relative")
      << AOpt2007(mVSharedIP,"SharedIP","Shared intrinsic parameters",{eTA2007::CanRepeat})



      << cHeaderSectionArg("Computation")
      << AOpt2007(mNbIter,"NbIter","Number of iterations",{eTA2007::HDV})
      << AOpt2007(mLVM,"LVM","Levenberg–Marquardt parameter (to have better conditioning of least squares)",{eTA2007::HDV})
      << AOpt2007(mViscPose,"PoseVisc","Sigma viscosity on pose")
      << AOpt2007(mShow_Cond,"Cond","Compute and show system condition number",{eTA2007::HDV})
      << AOpt2007(mParamShow_UK_UC,"UC_UK","Param for uncertainty & Show names of unknowns (tuning)")

      << AOpt2007(mCheckMeasureAdded,"CheckMeasureAdded","Do we check the adding of measures)",{eTA2007::Tuning,eTA2007::HDV})


           << cHeaderSectionArg("Blocks")
      << AOpt2007
         (
             mParamBOI,
             "BOI",
             "Bloc of Instr : [Pair-params],[Gauge-params],[Cur-rattach]?"
          )
      << AOpt2007
         (
             mParamBOIClino,
             "ClinoBOI",
             "Clino parameter"
          )
      << mPhProj.DPBlockInstr().ArgDirInOpt()
      << mPhProj.DPBlockInstr().ArgDirOutOpt()


    ;
}



void  cAppliBundlAdj::AddOneSetGCP3D(const std::string & aFolderIn, const std::string & aFolderOut, tREAL8 aWFactor, bool aDoExportSigmas)
{
    cSetMesGndPt  aFullMesGCP;
    cMes3DDirInfo * aMesDirInfo = cMes3DDirInfo::addMes3DDirInfo(mBA.getGCP(), aFolderIn, aFolderOut, aWFactor, aDoExportSigmas);
    mPhProj.LoadGCP3DFromFolder(aFolderIn, aFullMesGCP, aMesDirInfo, "", mGCPFilter, mGCPFilterAdd);
    auto aFullMes3D = aFullMesGCP.ExtractSetGCP("???");
    mBA.AddGCP3D(aMesDirInfo,aFullMes3D);
}


// VParam standar is done from  Folder +  weight of size [1,4]
void  cAppliBundlAdj::AddOneSetGCP2D(const cGCP2D & aGCP2D)
{
    mMeasureAdded = true;
    std::string aFolderIn = aGCP2D.GCPDir;  // folder
    std::vector<tREAL8>  aGCPW = {aGCP2D.Sigma, aGCP2D.SigmaAttenuation,aGCP2D.Threshold, aGCP2D.Exponent};
    cMes2DDirInfo * aMesDirInfo = cMes2DDirInfo::addMes2DDirInfo(mBA.getGCP() ,aFolderIn, cStdWeighterResidual(aGCPW,0));
    for (const auto  & aSens : mBA.VSIm())
    {
        // Load the images measure + init sens
        mPhProj.LoadImFromFolder(aFolderIn,mBA.getGCP().getMesGCP(),aMesDirInfo,aSens->NameImage(),aSens,SVP::Yes);
    }
}


void  cAppliBundlAdj::AddOneSetTieP(const std::string & aFolder, double aSigma, double aSigmaAttenuation, double aThreshold, double aExponent)
{
    mMeasureAdded = true;
    std::vector<tREAL8>  aTiePW = {aSigma, aSigmaAttenuation, aThreshold, aExponent};
    cStdWeighterResidual aWeighter(aTiePW,0);
    mBA.AddMTieP(aFolder,AllocStdFromMTPFromFolder(aFolder,VectMainSet(0),mPhProj,false,true,false),aWeighter);
}


int cAppliBundlAdj::Exe()
{

   // bool hasOrientPC     = false;
    bool hasConstrOriPC  = mLVM > 0;
/*
{
    StdOut() << "TESTT  mAddGCPWmAddGCPW\n";
    StdOut() << mAddGCPW  << mAddGCPW.size() << "\n";
    for (const auto& aGCP : mAddGCPW)
         StdOut() << " * " << aGCP  << aGCP.size() << "\n";
        


    getchar();
}
*/

    //   ========== [0]   initialisation of def values  =============================
 //  mPhProj.DPMulTieP().SetDirInIfNoInit(mDataDir);
  //  mPhProj.DPRigBloc().SetDirInIfNoInit(mDataDir); //  RIGIDBLOC

    mPhProj.FinishInit();

    {
        std::string aReportDir = mPhProj.DPOrient().DirIn() + "_" + mPhProj.DPOrient().DirOut();
        if (IsInit(& mPostFixReport))
           aReportDir += "_" + mPostFixReport;
        SetReportSubDir(aReportDir);
    }


    if (IsInit(&mParamRefOri))
    {
        hasConstrOriPC = true;
        mBA.AddReferencePoses(mParamRefOri.Ori, mParamRefOri.SigmaTr, mParamRefOri.SigmaRot, mParamRefOri.PatApply);
    }

    //   ========== [1]   Read unkowns of bundle  =============================
    for (const auto &  aNameIm : VectMainSet(0))
    {
         mBA.AddCam(aNameIm);
    }



    if (IsInit(&mPatParamFrozCalib))
    {
        mBA.SetParamFrozenCalib(mPatParamFrozCalib);
    }
    if (IsInit(&mVVParFreeCalib))
    {
        std::vector<cFreeCalibPattern> aVPat;
        for (const auto & aParam : mVVParFreeCalib)
            aVPat.push_back({aParam.PatCal, aParam.PatParam, aParam.Weight});
        mBA.SetParamFreeCalib(aVPat);
    }

    if (IsInit(&mPatFrosenCenters))
    {
        hasConstrOriPC = true;
        mBA.SetFrozenCenters(mPatFrosenCenters);
    }

    if (IsInit(&mPatFrosenOrient))
    {
        hasConstrOriPC = true;
        mBA.SetFrozenOrients(mPatFrosenOrient);
    }

    if (IsInit(&mViscPose))
    {
        hasConstrOriPC = true;
        mBA.SetViscosity(mViscPose.SigmaCenter,mViscPose.SigmaRot);
    }

    if (IsInit(&mVSharedIP))
    {
        mBA.SetSharedIntrinsicParams(mVSharedIP);
    }

    for (const auto& aGCP3D : mGCP3D)
    {
        hasConstrOriPC = true;

        AddOneSetGCP3D(aGCP3D.GCPDir,
                       aGCP3D.OutDir,
                       aGCP3D.Sigma,
                       aGCP3D.ExportSigma
                      );
    }

    if (mPhProj.DPTopoMes().DirInIsInit())
    {
        mMeasureAdded = true;
        mBA.AddTopo();
    }

    for (const auto& aGCP2D : mGCP2D)
    {
        // expected: [Folder,SigI,SigAt?=-1,Thrs?=-1,Exp?=1]
        AddOneSetGCP2D(aGCP2D);
    }

    if (mGCP2D.empty())
    {
        //here no 2d mes, fake it
        cMes2DDirInfo * aMes2DDirInfo = cMes2DDirInfo::addMes2DDirInfo(mBA.getGCP(), "in",cStdWeighterResidual());
        cSetMesPtOf1Im aSetMesPtOf1Im;
        mBA.AddGCP2D(aMes2DDirInfo, aSetMesPtOf1Im, nullptr, eLevelCheck::NoCheck);
    }


    if (IsInit(&mTiePWeight))
    {
        AddOneSetTieP(mPhProj.DPMulTieP().DirIn(), mTiePWeight.Sigma, mTiePWeight.SigmaAttenuation,
                      mTiePWeight.Threshold, mTiePWeight.Exponent);
    }

    mBA.SetTiePShowPerMil(mTiepShowPerMil);
    // Add  the potential suplementary TieP
    for (const auto& aTieP : mAddTieP)
        AddOneSetTieP(aTieP.Folder, aTieP.Sigma, aTieP.SigmaAttenuation, aTieP.Threshold, aTieP.Exponent);

    bool hasGauge = IsInit(&mParamGaugeRel);
    bool forceNoGauge =   (mParamGaugeRel.MainIm==MMVII_NONE);
    if ((!hasConstrOriPC) && (!hasGauge) && (mBA.NbCamPC()!=0) && (!forceNoGauge))
    {
            MMVII_USER_TYPED_WARNING(eTyUEr::eForceGauge,"Gauge in pure relative pause not specified, added by system");
            hasGauge = true;
    }

    if (hasGauge && (!forceNoGauge))
    {
        mBA.SetGaugeRelPause(mParamGaugeRel);
    }


    if (IsInit(&mParamLine))
    {
        mBA.AddLineAdjust(mParamLine.Folder, mParamLine.SigmaIm, mParamLine.NbPtsSampl);
    }

    if (IsInit(&mParamBOI))
    {
       mBA.AddBlockInstr(mParamBOI);
    }

    if (IsInit(&mParamBOIClino))
    {
       mBA.AddClinoBlokcInstr(mParamBOIClino);
    }
    


    if (IsInit(&mPatFrosenClino))
    {
        mBA.SetFrozenClinos(mPatFrosenClino);
    }
    
    for (const auto & aParam : mParamLidarPhgr)
    {
        mMeasureAdded = true;
        mBA.Add1AdjLidarPhotogra(aParam.Mode, aParam.PlyFile, aParam.Sigma, aParam.Interp,
                                 aParam.Pertubate, aParam.NbPtsPerPatch);
    }

    for (const auto & aParam : mParamLidarPhoto)
    {
        mMeasureAdded = true;
        mBA.Add1AdjLidarPhoto(aParam.Mode, aParam.PatScan, aParam.Sigma, aParam.Interp,
                              aParam.ScaleInit, aParam.ScaleFinal, aParam.Threshold, aParam.NbPtsPerPatch);
    }

    for (const auto & aParam : mParamLidarLidar)
    {
        mMeasureAdded = true;
        mBA.Add1AdjLidarLidar(aParam.PatScan, aParam.Sigma, aParam.ThresholdInit,
                              aParam.ThresholdFinal, aParam.NormalTolDeg, aParam.Interp);
    }

    if (mCheckMeasureAdded)
    {
        MMVII_INTERNAL_ASSERT_User(mMeasureAdded,eTyUEr::eUnClassedError,"Not any measure added");
    }

    if (IsInit(&mParamShow_UK_UC))
       mBA.Set_UC_UK(mParamShow_UK_UC);

    //   ========== [2]   Make Iteration =============================
    StdOut() << "Begin iterations" << std::endl;
    mBA.Iterate(mNbIter, mLVM, mShow_Cond);

    //   ========== [3]   Save resulst =============================
    for (auto & aSI : mBA.VSIm())
        mPhProj.SaveSensor(*aSI);
            /*
    for (auto & aCamPC : mBA.VSCPC())
        mPhProj.SaveCamPC(*aCamPC);
        */

    mPhProj.CpSysCoIn2Out(true,true);

    mBA.Save_newGCP3D();
    mBA.SaveTopo(); // just for debug for now
  //   mBA.SaveClino();
    mBA.SaveBlockInstr();

    if (IsInit(&mParamShow_UK_UC))
    {
        mBA.ShowUKNames(mParamShow_UK_UC,mPostFixReport,this);
    }

    return EXIT_SUCCESS;
}


tMMVII_UnikPApli Alloc_BundlAdj(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppliBundlAdj(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_OriBundlAdj
(
     "OriBundleAdj",
      Alloc_BundlAdj,
      "Bundle adjusment between images, using several observations/constraint",
      {eApF::Ori},
      {eApDT::Orient},
      {eApDT::Orient},
      __FILE__
);

/*
*/

}; // MMVII

