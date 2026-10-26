#include "MMVII_DeclareAllCmd.h"
#include "MMVII_DeclareCste.h"
#include "MMVII_Image2D.h"
#include "MMVII_MeasuresIm.h"
#include "MMVII_Random.h"
#include "MMVII_Sensor.h"

#include "../CodedTarget/CodedTarget.h"

#include <filesystem>


namespace MMVII
{


namespace fs = std::filesystem;
class ScopedChdir
{
    fs::path old_;
public:
    explicit ScopedChdir(const fs::path& dir) : old_(fs::current_path()) {
        fs::current_path(dir);
    }
    ~ScopedChdir() {
        std::error_code ec;
        fs::current_path(old_, ec);          // pas d'exception dans un destructeur
    }
    ScopedChdir(const ScopedChdir&) = delete;
    ScopedChdir& operator=(const ScopedChdir&) = delete;
};



class cPairTargetMatch
{
   public :
       cPairTargetMatch();

       std::vector<bool>    mDone;
       std::vector<cPt2dr>  mPts;
       tHom2Dr              mRefHT2Im;
       tAff2Dr              mExtAI2T;

};

cPairTargetMatch::cPairTargetMatch() :
    mDone ({false,false}),
    mPts  (2)
{

}


class cBenchCodedTarget
{
    public :
      cBenchCodedTarget(bool doOnlyVerif,bool showCom);
      void DoAll();

    private :

      /// Do all the tesst for a given
      void DoOneTest(eTyCodeTarget aType,int aNBB,int aNbTartget);

      /// [1] Generate random background image :
      void GenerateRandomImages(int aNumberOfImages,const cPt2di&  aSizeOfImages);

      /// [2] Generate the current file of specif
      void GenerateXmlTargetSpecificationFiles();
      /// [3]  Generate simulation
      void GenerateImagesWithTargets(const std::string & aNameFullSpec,int aNbTarget);
      /// [4] Run the target extraction programm adapted to current target model
      void DoTargetsExtraction(bool isIGN);

      /// [5]  do the verification for current target model
      void VerifTargets();
      void VerifTargetOf1Image(const std::string &);
      void  ReadMeasureIm(const std::string & aPrefix,const std::string & aNameIm,bool isRef);

      std::string NameRandomImage(int aKIm) const {return  mPrefixIm+ToStr(aKIm) + ".tif";}

      bool                                    mDoOnlyVerif;
      bool                                    mShowCom;
      int                                     mCurNbb;
      eTyCodeTarget                           mCurType;
      bool                                    mCurIsIgn;
      const cSpecMMVII_Appli *                mCurSpec;

      // ------------  different name used ------------------
      std::string                            mNameFullSpecif;
      const std::string                      mPrefixIm ;
      const std::string                      mPatImInit ;
      std::string                            mSuplPrefSim ;
      std::string                            mPatImSimul;


      //  -------------- Parameters for storing resiuduals -----------------
      std::map<std::string,cPairTargetMatch> mMapPairs;
      cStdStatRes                            mPtsRes;
      int                                    mNbMissedDet;
      int                                    mNbOverDet;
};

cBenchCodedTarget::cBenchCodedTarget(bool doOnlyVerif,bool showCom) :
     mDoOnlyVerif (doOnlyVerif),
     mShowCom     (showCom),
     mPrefixIm    ( "ImSynt_"),
     mPatImInit   (mPrefixIm+".*tif")
{
}




/**
 * @brief cBenchCodedTarget::GenerateImage
 * @param aNbIm
 * @param aSzIm
 *
 *   Generate  NbIm background image for target incrustation
 */
void   cBenchCodedTarget::GenerateRandomImages(int aNbIm,const cPt2di&  aSzIm)
{
    for (int aKIm=0 ; aKIm<aNbIm ; aKIm++)
    {
        // MMVII ImageGenRandom "[1500,1000]" Type=U_INT1 Out=ImSynt_0.tif
        //    SeedRand=0 "GaussNoise=[100,70]" "GaussNoise=[5,10]" StdOut=NONE

        cParamCallSys aCom(TheSpec_cAppliGenRandomImage,true);
        aCom.AddArgs(ToStr(aSzIm));

        aCom.AddMMVIIArgsOpt("Type","U_INT1");  // Make 8-bits
        aCom.AddMMVIIArgsOpt(CurOP_Out,NameRandomImage(aKIm));  // Generate a diff image for each K

        // each image is different , set the seed, +- majic formula using prime number , Chinese Remainder Theorem
        int aSeed = 11 + cRandGenerator::TheOne()->GetCurrentSeed() * 1009 + 37*aKIm;
        aCom.AddMMVIIArgsOpt("SeedRand",ToStr(aSeed));

        std::vector<std::vector<double>> aVecNoise{{100.0,70.0},{5.0,10.0}};
        for (const auto & aGaussN : aVecNoise)
            aCom.AddMMVIIArgsOpt("GaussNoise",ToStr(aGaussN));

        aCom.ExecuteWithoutMsg(mShowCom);
     }
}


void  cBenchCodedTarget::GenerateXmlTargetSpecificationFiles()
{

    // Opporunity to test dmp + faster when many target
    std::string aPost = (mCurNbb<=14 ) ? PostF_XmlFiles : PostF_DumpFiles ;

    //  ----------------- Generate the encoding -----------------------------------------------

      std::string aNameEncoding = "Encoding_" + E2Str(mCurType) + "_" + ToStr(mCurNbb) + "." + aPost;
      {
          //  MMVII CodedTargetGenerateEncoding IGNIndoor 14 Out=Encoding_IGNIndoor_14.xml StdOut=NONE

          cParamCallSys aComEncoding(TheSpecGenerateEncoding,true);
          aComEncoding.AddArgs(E2Str(mCurType));
          aComEncoding.AddArgs(ToStr(mCurNbb));
          aComEncoding.AddMMVIIArgsOpt(CurOP_Out,aNameEncoding);

          aComEncoding.ExecuteWithoutMsg(mShowCom);
      }

      //  ----------------- Generate the full spec, with geometry  -----------------------------------------------

      mNameFullSpecif = "FullSpec_" + E2Str(mCurType) + "_" + ToStr(mCurNbb) + "." + aPost;
      {
          // MMVII CodedTargetGenerate Encoding_IGNIndoor_14.xml Out=FullSpec_IGNIndoor_14.xml StdOut=NONE

          cParamCallSys aComGenerate(TheSpecGenCodedTarget,true);
          aComGenerate.AddArgs(aNameEncoding);
          aComGenerate.AddMMVIIArgsOpt(CurOP_Out,mNameFullSpecif);

          // For (still) unexplained reason, there is a memory leak if we do  Execute(false)
          // (do it in the same process) by the way, the command if it goes untill the end is probably
          // memory correct as its has been widely used . So, for now, run in a separate process

          aComGenerate.ExecuteWithoutMsg(mShowCom);
      }
}


void cBenchCodedTarget::GenerateImagesWithTargets(const std::string & ,int aNbTarget)
{

    // MMVII CodedTargetSimul "ImSynt_.*tif" FullSpec_IGNIndoor_14.xml SuplPref=IGNIndoor14_ NbSel=15
    //  "Radius=[40,60]" "Ratio=[0.7,1]" "NoiseAmpl=[0,0]" "PropLinBias=[0,0]" StdOut=NONE

    cParamCallSys aComSimul(TheSpecSimulCodedTarget,true);

    aComSimul.AddArgs(mPatImInit);
    aComSimul.AddArgs(mNameFullSpecif);

    aComSimul.AddMMVIIArgsOpt("SuplPref",mSuplPrefSim);
    aComSimul.AddMMVIIArgsOpt("NbSel",ToStr(aNbTarget));
    aComSimul.AddMMVIIArgsOpt("Radius",ToStr(cPt2dr(40,60)));
    aComSimul.AddMMVIIArgsOpt("Ratio",ToStr(cPt2dr(0.7,1.0)));
    aComSimul.AddMMVIIArgsOpt("NoiseAmpl",ToStr(cPt2dr(0,0)));
    aComSimul.AddMMVIIArgsOpt("PropLinBias",ToStr(cPt2dr(0,0)));


    aComSimul.ExecuteWithoutMsg(mShowCom);
}


void cBenchCodedTarget::DoTargetsExtraction(bool )
{
    // For ex with IGN target :
    //
    //    MMVII CodedTargetCheckBoardExtract "SimulTarget_IGNIndoor14_ImSynt_.*tif"
    //        FullSpec_IGNIndoor_14.xml OutObjMesInstr=Bench StdOut=NONE

    cParamCallSys aComExtract(*mCurSpec,true);

    aComExtract.AddArgs(mPatImSimul);
    aComExtract.AddArgs(mNameFullSpecif);

    aComExtract.AddMMVIIArgsOpt(eTA2007::ObjMesInstr,false,"Bench");
   // aComExtract.AddMMVIIArgsOpt("OutObjMesInstr","Bench");

    aComExtract.ExecuteWithoutMsg(mShowCom);
}


//cSetMesPtOf1Im  FileMeasureIm(const std::string & aPrefix,const std::string & aNameIm)
void  cBenchCodedTarget::ReadMeasureIm(const std::string & aPrefix,const std::string & aNameIm,bool isRef)

{
  /*  std::string aDir =     MMVII_DirPhp
                         + E2Str(eTA2007::ObjMesInstr) +  StringDirSeparator()
                         + aPrefix +   StringDirSeparator();*/

    std::string aDir = cDirsPhProj::StdDirOfType(eTA2007::ObjMesInstr,aPrefix);
    std::string aRes = aDir + cSetMesPtOf1Im::StdNameFileOfIm(aNameIm);;



    {
        cSetMesPtOf1Im aSet =  cSetMesPtOf1Im::FromFile(aRes);

        for (const auto &  aMes : aSet.Measures())
        {
            cPairTargetMatch & aPair =  mMapPairs[aMes.mNamePt];
            MMVII_INTERNAL_ASSERT_always(! aPair.mDone[isRef],"Unexpected case in cBenchCodedTarget::ReadMeasureIm");
            aPair.mDone[isRef] = true;
            aPair.mPts[isRef] = aMes.mPt;
        }
    }


    if (isRef)
    {
        std::string aName = aDir + NameGTSimulTargetFromIm(aNameIm);
        cResSimul aRS;
        ReadFromFile(aRS,aName);

        for (const auto & aGT : aRS.mVG)
        {
            std::string aName = aGT.mEncod.Name();
            cPairTargetMatch & aPair =  mMapPairs[aName];
            MMVII_INTERNAL_ASSERT_always(aPair.mDone[isRef],"Incoh in cBenchCodedTarget::ReadMeasureIm");
            aPair.mRefHT2Im = aGT.mHomogT2Im;
        }

    }
    else
    {
        std::string aName = aDir + cSaveExtrEllipe::PrefixAttr + "-" +  cSetMesPtOf1Im::StdNameFileOfIm(aNameIm);

        std::vector<cSaveExtrEllipe> aVecSEE;
        ReadFromFile(aVecSEE,aName);
        for (const auto & aSEE : aVecSEE )
        {
            std::string aName = aSEE.mNameCode;
            cPairTargetMatch & aPair =  mMapPairs[aName];
            MMVII_INTERNAL_ASSERT_always(aPair.mDone[!isRef],"Incoh in cBenchCodedTarget::ReadMeasureIm");
            aPair.mExtAI2T = aSEE.mAffIm2Ref;
        }
        StdOut() << "MMM= "  << aName << " Sz="<< aVecSEE.size() <<  "\n";

    }

}



void cBenchCodedTarget::VerifTargetOf1Image(const std::string &aNameIm)
{
    mMapPairs.clear();
    ReadMeasureIm("Simul",aNameIm,true);
    ReadMeasureIm("Bench",aNameIm,false);

    for (const auto & [aName,aPair] : mMapPairs)
    {
        if ((!aPair.mDone[0]) && aPair.mDone[1])
        {
           mNbOverDet++;
        }
        else if (aPair.mDone[0] && (!aPair.mDone[1]) )
        {
            mNbMissedDet++;

        }
        else  if (aPair.mDone[0] && aPair.mDone[1])
        {
            tREAL8 aRes = Norm2(aPair.mPts[0]-aPair.mPts[1]);
            mPtsRes.Add(aRes);

            cPt2dr aC = aPair.mRefHT2Im.Inverse(aPair.mPts[1]);
            StdOut() << "DDD=" << Norm2(aPair.mExtAI2T.Inverse(aC) - aPair.mRefHT2Im.Value(aC)) << "\n";
        }
        else
        {
            MMVII_INTERNAL_ERROR("Unexpected case in cBenchCodedTarget::VerifTargetOf1Image");
        }
    }
}

void cBenchCodedTarget::VerifTargets()
{
    mNbMissedDet=0;
    mNbOverDet=0;
    mPtsRes=    cStdStatRes();


    std::vector<std::string> aVecIm =  GetFilesFromDir("./",AllocRegex(mPatImSimul));
    // StdOut() << "VECIMSIM=" << aVecIm << "\n";

    for (const auto & aNameIm : aVecIm)
       VerifTargetOf1Image(aNameIm);

    StdOut() << E2Str(mCurType) << "/" << mCurNbb << " :: "
             << " Missed=" << mNbMissedDet
             << " Over=" << mNbOverDet
             << " Res=" << mPtsRes.Avg()
             << "\n";
}

void  cBenchCodedTarget::DoOneTest(eTyCodeTarget aType,int aNBB,int aNbTartget)
{
    mCurType = aType;
    mCurIsIgn = (mCurType!= eTyCodeTarget::eCERN);
    mCurSpec=  (mCurIsIgn ? &TheSpecExtractCheckBoardTarget : &TheSpecExtractCircTarget);

    mCurNbb = aNBB;
    mSuplPrefSim = E2Str(aType)+ToStr(aNBB) +"_";
    mPatImSimul  =  ThePrefixSimulTarget + mSuplPrefSim+mPrefixIm +".*tif";

    if (! mDoOnlyVerif)
    {
        GenerateXmlTargetSpecificationFiles();
        GenerateImagesWithTargets(mNameFullSpecif,15);
        DoTargetsExtraction(aType!=eTyCodeTarget::eCERN);
    }

    VerifTargets();

}


void cBenchCodedTarget::DoAll()
{
    std::string aDir = cMMVII_Appli::TmpDirTestMMVII() ;

    if (UserIsMPD())
    {
        aDir = "/home/MPierrot-Deseilligny/MMVII/MMVII-TestDir/Tmp/";
        StdOut() << " DIR=" << aDir << "\n";
    }

    // push to the tmp directory
    ScopedChdir aScopeChD(aDir);


    if (! mDoOnlyVerif)
    {
        GenerateRandomImages(3,cPt2di(1500,1000));
    }

     DoOneTest(eTyCodeTarget::eIGNIndoor,14,15);
    // DoOneTest(eTyCodeTarget::eCERN,14,15);
    // DoOneTest(eTyCodeTarget::eCERN,20,15);


}

void BenchCodedTarget(cParamExeBench & aParam)
{


    if (! aParam.NewBench("CodedTarget"))
        return;

    if  (!UserIsMPD())
    {
        aParam.EndBench();
        return;
    }
    else
    {
        MMVII_DEV_WARNING("BenchCodedTarget not active w/o MPD");
    }
    bool OnlyVerif = AppliDevloperGetBoolSpecValue("OnlyVerif","If true, do only the comparison extract/ref ");
    bool ShowCom = AppliDevloperGetBoolSpecValue("ShowCom","If true, print the comment executed");

    if (ShowUnusedDevKey())
    {
        aParam.EndBench();
        return;
    }

    StdOut() << " SHOW COM=" << ShowCom << "\n";


    cBenchCodedTarget aBCT(OnlyVerif,ShowCom);

    aBCT.DoAll();

    aParam.EndBench();
}



};
