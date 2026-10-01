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
            std::pair<cParamCallSys,cParamCallSys> StrComMasqMaker(const tNamePair * aPair);
            cParamCallSys StrComDenseMatchEpip(const tNamePair * aPair, int aGpuId);
            cParamCallSys StrComProjGround(const tNamePair * aPair);
        private :
            cPhotogrammetricProject  mPhProj;
            std::string  mPatternIm;
            std::string mPairsFile;
            tNameRel mNamePairs; // List of pairs of images to match, if empty pairs will be computed given orientation and tie points
            std::vector<const tNamePair *> mVecPairs;
            std::string mDirSelectPairs;
            std::string mDirEpipolarResampling;
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

    std::pair<cParamCallSys,cParamCallSys> cAppli_DenseMatch::StrComMasqMaker(const tNamePair * aPair)
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
        
        return std::make_pair(
            cParamCallSys (
                MMV1Bin(),
                "MasqMaker",
                mDirEpipolarResampling + "/" + aNameEpipIm1,
                "0",
                "255"
            ),
            cParamCallSys (
                MMV1Bin(),
                "MasqMaker",
                mDirEpipolarResampling + "/" + aNameEpipIm2,
                "0",
                "255"
            )
        );
    }

    cParamCallSys cAppli_DenseMatch::StrComDenseMatchEpip(const tNamePair * aPair, int aGpuId)
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
            "DirMEC="+mDirDenseMatchingEpip+"/"+aNameDMPerPair+"/",
            "DoCorrel=1",
            "DirProj="+mDirEpipolarResampling,
            "NbProc=1" // 3 proc per GPU
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
            "GroundResolution="+ToStr(mGSD)
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
        std::list<cParamCallSys> aVecComMasqMaker;
        std::list<cParamCallSys> aVecComDenseMatch;
        std::list<cParamCallSys> aVecComProjGround;

        int aGpuId = 0;
        if (mSetGpuPool.size() == 0)
            aGpuId = -1; // No GPU

        int aDBG_MAX =5;
        
        for (const auto & aPair : mVecPairs)
        {
            cParamCallSys aComEpipResample = StrComEpipResample(aPair);
            std::pair<cParamCallSys,cParamCallSys> aComMasqMaker = StrComMasqMaker(aPair);
            
            cParamCallSys aComDenseMatch   = StrComDenseMatchEpip(aPair,aGpuId);
            cParamCallSys aComProjGround   = StrComProjGround(aPair);

            if (aDBG_MAX>0)
            {
                StdOut() <<  aComEpipResample.Com() << "\n";
                aDBG_MAX--;
                aVecComEpipResample.push_back(aComEpipResample);

                // make mask of defined pixels after epipolar resampling
                aVecComMasqMaker.push_back(aComMasqMaker.first);
                aVecComMasqMaker.push_back(aComMasqMaker.second);
                StdOut() <<  aComMasqMaker.first.Com() << "\n";
                StdOut() <<  aComMasqMaker.second.Com() << "\n";
                // dense matching vect
                aVecComDenseMatch.push_back(aComDenseMatch);
                StdOut() <<  aComDenseMatch.Com() << "\n";
                if (mSetGpuPool.size() > 0)
                    aGpuId = (aGpuId + 1) % mSetGpuPool.size();
                aVecComProjGround.push_back(aComProjGround);
                StdOut() <<  aComProjGround.Com() << "\n";
            }
            else
            {
                break;
            }
        }

        // run Epipolar resampling
        //ExeComParal(aVecComEpipResample,true);
        //run dense matching
        //ExeComSerial(aVecComDenseMatch,true);
        //project to ground
        //ExeComSerial(aVecComProjGround,true);

        //Finally, fusion of all clouds
        cParamCallSys aComCloudFusion(
            "MMVII",
            "CloudMMVII_Fuse",
            mDirCloudProjOnGround+"Prof.*xml",
            "DirProj="+mDirCloudProjOnGround
        );
        //ExeComSerial({aComCloudFusion},true);
        StdOut() <<  aComCloudFusion.Com() << "\n";
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



