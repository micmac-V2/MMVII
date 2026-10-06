#include "cMMVII_Appli.h"
#include "MMVII_PCSens.h"
#include "MMVII_AllClassDeclare.h"
#include "MMVII_Tpl_ElemStrToVal.h"
 


namespace MMVII
{

    class cAppli_DenseMatch : public cMMVII_Appli
    {
        public :
            cAppli_DenseMatch(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec);
            int Exe() override;
            cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
            cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;
            cParamCallSys StrComEpipResample(const tNamePair * aPair);
            cParamCallSys StrComDenseMatchEpip(const tNamePair * aPair, int aGpuId, int aNbProc);
            cParamCallSys StrComProjGround(const tNamePair * aPair);
        private :
            cPhotogrammetricProject  mPhProj;
            std::string  mPatternIm;
            std::string mPairsFile;
            tNameRel mNamePairs; // List of pairs of images to match, if empty pairs will be computed given orientation and tie points
            std::vector<const tNamePair *> mVecPairs;
            std::string mDirSelectPairs;
            std::string mDirEpipolarResampling;
            bool mExec;  // If 1, execute the command, if 0, only print it
            std::string mDirDenseMatchingEpip;
            std::string mDirCloudProjOnGround;
            std::string mDirCloudFusion;
            std::string mNameEpipPattern;
            std::vector<int> mSetGpuPool;
            tREAL8 mGSD;  
    };


    cAppli_DenseMatch::cAppli_DenseMatch(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec) :
        cMMVII_Appli  (aVArgs,aSpec),
        mPhProj       (*this),
        mPairsFile    (""),
        mDirSelectPairs(""),
        mDirEpipolarResampling(""),
        mDirDenseMatchingEpip(""),
        mDirCloudProjOnGround(""),
        mDirCloudFusion(""),
        mNameEpipPattern("Epip_%1_%2.tif"),
        mGSD(0.1)
    {
    }   

    cCollecSpecArg2007 & cAppli_DenseMatch::ArgObl(cCollecSpecArg2007 & anArgObl)
    {
        return anArgObl
               << Arg2007(mPatternIm,"Pattern for images",{{eTA2007::MPatFile,"0"}})
               << mPhProj.DPOrient().ArgDirInMand()
            ;
    }

    cCollecSpecArg2007 & cAppli_DenseMatch::ArgOpt(cCollecSpecArg2007 & anArgOpt)
    {
        return anArgOpt
            << AOpt2007(mPairsFile,"FilePairs","File containing pairs of images to match",{eTA2007::FileTxt})
            <<  mPhProj.DPTieP().ArgDirInOpt()
            <<  mPhProj.DPMulTieP().ArgDirInOpt()
            << AOpt2007(mSetGpuPool,"GpusList", "List of Gpu Ids to dispatch computation if ", {eTA2007::CanRepeat})
            << AOpt2007(mGSD,"GndSamplingD", "Ground sampling distance of bascule")
            << AOpt2007(mDirEpipolarResampling,"DirDenseMatch","Directory for overall Dense Matching")
            << AOpt2007(mExec,"Exec","If 1, execute the command, if 0, only print it")
            ;
    }


    cParamCallSys cAppli_DenseMatch::StrComEpipResample(const tNamePair * aPair)
    {
        
        std::string aNameIm1 = aPair->V1();
        std::string aNameIm2 = aPair->V2();

        return cParamCallSys (
            "MMVII",
            "EpipResampling",
            aNameIm1,
            aNameIm2,
            mPhProj.DPOrient().DirIn(),
            mPhProj.DPTieP().DirInIsInit() ? "TieP=" + mPhProj.DPTieP().DirIn() : "MulTieP=" + mPhProj.DPMulTieP().DirIn(),
            "TiePMinNbRatio="+ToStr(0.006)
            //"OutDir="+mPhProj.DirPhp()+"/VISU/DenseMatch/"+aNameIm1+"_"+aNameIm2
        );
    }

    cParamCallSys cAppli_DenseMatch::StrComDenseMatchEpip(const tNamePair * aPair, int aGpuId, int aNbProc)
    {
        std::string aName1 = LastPrefix(aPair->V1());
        std::string aName2 = LastPrefix(aPair->V2());

        std::string aNameEpipIm1 = replaceFirstOccurrence(
                            replaceFirstOccurrence(mNameEpipPattern,"%1",aName1),
                            "%2",
                            aName2);
        std::string aNameEpipIm2 = replaceFirstOccurrence(
                            replaceFirstOccurrence(mNameEpipPattern,"%1",aName2),
                            "%2",
                            aName1);
        
        std::string aNameDMPerPair = "RAFTStereo_" +aNameEpipIm1;

        return cParamCallSys (
            "MMVII",
            "DenseMatchEpipGen",
            "RAFTStereo",
            aNameEpipIm1,
            aNameEpipIm2,
            "OnGPU="+ToStr(aGpuId),
            "DirMEC="+aNameDMPerPair+"/",
            "DoCorrel=1",
            "DirProj="+mDirEpipolarResampling,
            "NbProc="+ToStr(aNbProc)
        );
    }

    cParamCallSys cAppli_DenseMatch::StrComProjGround(const tNamePair * aPair)
    {
        std::string aName1 = LastPrefix(aPair->V1());
        std::string aName2 = LastPrefix(aPair->V2());

        std::string aNameEpipIm1 = replaceFirstOccurrence(
                            replaceFirstOccurrence(mNameEpipPattern,"%1",aName1),
                            "%2",
                            aName2);

        std::string aNameEpipIm1Masq = LastPrefix(aNameEpipIm1)+"_Masq.tif";

        std::string aNameEpipIm2 = replaceFirstOccurrence(
                            replaceFirstOccurrence(mNameEpipPattern,"%1",aName2),
                            "%2",
                            aName1);
        
        std::string aNameDMPerPair = "RAFTStereo_" +aNameEpipIm1;

        return cParamCallSys (
            "MMVII",
            "CloudProjectOnGround",
            aNameEpipIm1,
            aNameDMPerPair+"/Px1_Num3_DeZoom1_LeChantier.tif",
            "Epi",
            "BASC",
            "Masq1="+aNameEpipIm1Masq,
            "Im2="+aNameEpipIm2,
            "ImCorrel="+aNameDMPerPair+"/Correl_LeChantier_Num3.tif",
            "DirProj="+mDirCloudProjOnGround,
            "GroundResolution="+ToStr(mGSD),
            "OriInFatherDir=1",
            "NbProc=1"
        );
    }

    int cAppli_DenseMatch::Exe()
    {
        mPhProj.FinishInit();

        mDirSelectPairs= mPhProj.DirVisu()+ "/" + "DMSelectBestPairs"+ "/";
        mDirEpipolarResampling=mPhProj.DirVisu() + "/" + "EpipResampling" + "/";
        mDirDenseMatchingEpip=mDirEpipolarResampling;
        mDirCloudProjOnGround=mDirEpipolarResampling;
        mDirCloudFusion=mDirEpipolarResampling;

        StdOut() << "DirVisu=" << mPhProj.DirVisu() << "\n";
        //1. Compute Image pairs if not given
        if (IsInit(&mPairsFile))
        {
            // read pairs file
            bool IsExit = ExistFile(mPairsFile);
            MMVII_INTERNAL_ASSERT_User(IsExit, eTyUEr::eUnClassedError , "Pairs file does not exist");

            mNamePairs = RelNameFromXmlFileIfExist(mPairsFile,IsExit)   ;

            StdOut()  << "Read " << mNamePairs.size() << " pairs of images from file: " << mPairsFile << "\n";  

            mNamePairs.PutInVect(mVecPairs,true);  // Sort the pairs, to have a deterministic order
            MMVII_INTERNAL_ASSERT_User((mNamePairs.size()>0), eTyUEr::eUnClassedError , "No pairs of images to match");
        }
        else
        {
            // Compute with MMVII DMSelectBestPairs
           std::string aWhichTieP= mPhProj.DPTieP().DirInIsInit() ? "InTieP=" + mPhProj.DPTieP().DirIn() : "InMulTieP=" + mPhProj.DPMulTieP().DirIn();
           cParamCallSys aComComputePairs(
                "MMVII",
                "DMSelectBestPairs",
                mPatternIm,
                mPhProj.DPOrient().DirIn(),
                aWhichTieP,
                "NbMinHomol=10",
                "NbMaxHomol=1000",
                "CellSize=10"
          );
          
           aComComputePairs.Execute(true);

           bool IsExit = ExistFile(mDirSelectPairs+"/"+"Pairs.xml");

           mNamePairs = RelNameFromXmlFileIfExist(mDirSelectPairs+"/"+"Pairs.xml",IsExit)   ;

           mNamePairs.PutInVect(mVecPairs,true);  // Sort the pairs, to have a deterministic order
           
           MMVII_INTERNAL_ASSERT_User((mNamePairs.size()>0), eTyUEr::eUnClassedError , "No pairs of images to match");

        }

        //2. For each pair,  epipolar rectificaton, dense matching, projection to ground
        std::list<cParamCallSys> aVecComEpipResample;
        std::list<cParamCallSys> aVecComDenseMatch;
        std::list<cParamCallSys> aVecComProjGround;

        int aGpuId = 0;
        int aNbProc =3 ;
        if (mSetGpuPool.size() == 0)
            aGpuId = -1; // No GPU
        
        for (const auto & aPair : mVecPairs)
        {
            cParamCallSys aComEpipResample = StrComEpipResample(aPair);
            aNbProc = ( aGpuId==-1 ) ? 1 : aNbProc;
            cParamCallSys aComDenseMatch   = StrComDenseMatchEpip(aPair,aGpuId,aNbProc);
            cParamCallSys aComProjGround   = StrComProjGround(aPair);

            aVecComEpipResample.push_back(aComEpipResample);
            StdOut()<< aComEpipResample.Com()<<std::endl;
            // dense matching vect
            aVecComDenseMatch.push_back(aComDenseMatch);
            StdOut() <<  aComDenseMatch.Com() << "\n";
            if (mSetGpuPool.size() > 0)
                aGpuId = (aGpuId + 1) % mSetGpuPool.size();
            aVecComProjGround.push_back(aComProjGround);
            StdOut() <<  aComProjGround.Com() << "\n";
        }

        //Finally, fusion of all clouds
        cParamCallSys aComCloudFusion(
            "MMVII",
            "CloudMMVII_Fuse",
            "Prof.*xml",
            "DirProj="+mDirCloudProjOnGround
        );

        StdOut() <<  aComCloudFusion.Com() << "\n";

        if (mExec)
        {
            // run Epipolar resampling
            ExeComParal(aVecComEpipResample,true);
            //run dense matching
            // set allowed proc to 3
            mNbProcAllowed = 3;
            ExeComSerial(aVecComDenseMatch,true);
            //project to ground
            // large nbproc for projonground 
            mNbProcAllowed=12;
            // later scale with computation load
            ExeComParal(aVecComProjGround,true);
            ExeComSerial({aComCloudFusion},true);
        }
        else
        {
             StdOut() << "Execution disabled, only printing commands\n";
        }

        return EXIT_SUCCESS;
    }

    /*  ============================================= */
    /*       ALLOCATION                               */
    /*  ============================================= */

    tMMVII_UnikPApli Alloc_DenseMatch(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec)
    {
        return tMMVII_UnikPApli(new cAppli_DenseMatch(aVArgs,aSpec));
    }

    cSpecMMVII_Appli  TheSpecDenseMatch
        (
            "DenseMatch",
            Alloc_DenseMatch,
            "Dense matching of a set of images, given orientation and tie points",
            {eApF::Match},
            {eApDT::Image,eApDT::Ori},
            {eApDT::Image},
            __FILE__
            );

};



