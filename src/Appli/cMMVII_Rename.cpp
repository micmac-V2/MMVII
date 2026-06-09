#include "cMMVII_Appli.h"
#include "MMVII_Sys.h"
#include "MMVII_DeclareCste.h"
#include "MMVII_2Include_Serial_Tpl.h"
#include "MMVII_Sensor.h"

#include <regex>
#include <filesystem>


namespace MMVII
{
std::string ArithmReplace(const std::string & aStrIn0,const std::regex& aPat,const std::vector<std::string> & aVArg)
{
    std::string anOp = aVArg[0];
    int aOffset = cStrIO<int>::FromStr(aVArg[1]);
    int aKExpr = cStrIO<int>::FromStr(aVArg[2]);

    std::smatch aBoundMatch;
    bool aGotMatch = std::regex_search(aStrIn0, aBoundMatch, aPat);
    MMVII_INTERNAL_ASSERT_tiny(aGotMatch,"cCRegex::BoundsMatch no match");

  //  StdOut() << " KKK=" << aKExpr  << " BM=" << aBoundMatch.size() << "\n";

    if ((aKExpr<0)||(aKExpr >= (int) aBoundMatch.size()))
    {
     //   StdOut() << " PAT=" << <<"\n";
         MMVII_UnclasseUsEr("Num of expr incompatible with pattern : " );
    }

    // auto aMatch  = aBoundMatch[aKExpr];

    std::string aStrNumIn = aBoundMatch[aKExpr];

  //  StdOut() << " NNNN [" << aStrNumIn << "]\n";
    int aNum    = cStrIO<int>::FromStr(aStrNumIn);

    if (anOp=="+")
       aNum += aOffset;
    else if (anOp=="-")
       aNum -= aOffset;
    else if (anOp=="%")
       aNum %= aOffset;
    else
    {
       MMVII_UnclasseUsEr("Bad operand in arithmetic : " + aVArg[0]);
    }

    int aNbDig = (aVArg.size()> 3) ? cStrIO<int>::FromStr(aVArg[3]) : aStrNumIn.size();
    std::string aStrNumOut = ToStr(aNum,aNbDig);

    std::string aStrIn = aStrIn0;
    aStrIn.replace(aBoundMatch.position(aKExpr),aBoundMatch.length(aKExpr),aStrNumOut);

    return aStrIn;
}

/*
void CreateLink(const std::string & aFileTarget,const std::string & aLink2Create,bool fileMustExist = true)
{

    if (std::filesystem::is_symlink(aLink2Create))
    {
        std::string aPrevTarget  = std::filesystem::read_symlink(aLink2Create).native();
        if (aPrevTarget == aFileTarget)
        {
            MMVII_USER_WARNING("Link already exist pointing to same file :" + aLink2Create + "->" + aFileTarget );
        }
        else
        {
             MMVII_USER_WARNING
             (
                  "Link already exist pointing to diff file, do noting remove before :"
                + aLink2Create + "->" +  aPrevTarget + "/" + aFileTarget
             );
        }
        return;
    }
    if (fileMustExist)
    {
        MMVII_INTERNAL_ASSERT_always(ExistFile(aFileTarget),"File "+aFileTarget + " dont exist in CreateLink");
    }
    std::filesystem::create_symlink(aFileTarget,aLink2Create);
}
*/

// std::filesystem::path read_symlink( const std::filesystem::path& p );
// bool is_symlink( const std::filesystem::path& p );

/* ==================================================== */
/*                                                      */
/*          cAppli_DicoRename                           */
/*                                                      */
/* ==================================================== */


class cAppli_DicoRename   : public cMMVII_Appli
{
     public :
        cAppli_DicoRename  (const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli &);  ///< constructor
        int Exe() override;                                             ///< execute action
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override; ///< return spec of  mandatory args
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override; ///< return spec of optional args
     protected :
     private :
        cPhotogrammetricProject  mPhProj;
        std::string              mNameFileTxtIn;
        std::string              mFormat;
        std::vector<std::string> mPatIm;
        std::string              mNameDico;
        int                      mL0;
        int                      mLLast;
        char                     mComment;
        std::string              mSeparator;
        int                      mNbMinTieP;
        std::vector<std::string> mNameFilesListIm;

};

cAppli_DicoRename::cAppli_DicoRename(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
    cMMVII_Appli      (aVArgs,aSpec),
    mPhProj           (*this),
    mL0               (0),
    mLLast            (-1),
    mComment          ('#'),
    mSeparator        ("@"),
    mNbMinTieP        (0),
    mNameFilesListIm  {"AllImDicoIn.xml","AllImDicoOut.xml"}
{
}

cCollecSpecArg2007 & cAppli_DicoRename::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return anArgObl
              <<  Arg2007(mNameFileTxtIn ,"Name of Input File")
              <<  Arg2007(mFormat   ,"Format of file as for ex \"SNSXYZSS\" ")
              <<  Arg2007(mPatIm ,"Substitution pattern [Pattern,SubstIn,SubstOut]",{{eTA2007::ISizeV,"[3,3]"}})
              <<  Arg2007(mNameDico ,"Name for output dictionnary")
           ;
}

cCollecSpecArg2007 & cAppli_DicoRename::ArgOpt(cCollecSpecArg2007 & anArgObl)
{

    return anArgObl
       << AOpt2007(mL0,"NumL0","Num of first line to read",{eTA2007::HDV})
       << AOpt2007(mLLast,"NumLast","Num of last line to read (-1 if at end of file)",{eTA2007::HDV})
       << AOpt2007(mComment,"Com","Carac for commentary",{eTA2007::HDV})
       << AOpt2007(mNameFilesListIm,"Files","Name file to transform [Input,Output]",{{eTA2007::ISizeV,"[2,2]"},eTA2007::HDV})
       << AOpt2007(mNbMinTieP,"NbMinTiep","Number minimal of tie point for save, set -1 if save w/o tiep",{eTA2007::HDV})


       <<   mPhProj.DPMulTieP().ArgDirInOpt()
       <<   mPhProj.DPMulTieP().ArgDirOutOpt()

       <<   mPhProj.DPGndPt2D().ArgDirInOpt()
       <<   mPhProj.DPGndPt2D().ArgDirOutOpt()
     ;
}


int cAppli_DicoRename::Exe()
{
    mPhProj.FinishInit();

    std::vector<std::vector<std::string>> aVVNames;
    std::vector<std::vector<double>> aVNums;
    std::vector<cPt3dr> aVXYZ,aVWKP;

    ReadFilesStruct
    (
        mNameFileTxtIn, mFormat,
        mL0, mLLast, mComment,
        aVVNames,aVXYZ,aVWKP,aVNums,
        false
    );

    std::map<std::string,std::string>  aDico;
    for (auto & aVNames  : aVVNames)
    {
         std::string aCatName = aVNames.at(0);
         for (size_t aKName=1 ; aKName<aVNames.size() ; aKName++)
             aCatName = aCatName + mSeparator + aVNames.at(aKName);

         std::string  aNameIn  = ReplacePattern(mPatIm.at(0),mPatIm.at(1),aCatName);
         std::string  aNameOut = ReplacePattern(mPatIm.at(0),mPatIm.at(2),aCatName);

         aDico[aNameIn] = aNameOut;

         // StdOut()  << "DddDgkI :: " << aNameIn  << " => " << aNameOut << "\n";
    }

    SaveInFile(aDico,mNameDico);

/*
    if (IsInit(&mNameFiles))
    {
        auto aSetIn = ToVect(SetNameFromString(mNameFiles.at(0),true));
        tNameSet aSetIn;
        tNameSet aSetOut;

        for (const auto & aNameIn : aSetIn)
        {
            const auto & anIter = aDico.find(aNameIn);
            if (anIter != aDico.end())
            {
               aSetOut.Add(anIter->second);
            }
        }
        SaveInFile(aSetIn,mNameFiles.at(0));
        SaveInFile(aSetOut,mNameFiles.at(1));
    }
*/

    bool  isInitMTP = mPhProj.DPMulTieP().DirInIsInit();
    if (isInitMTP)
    {
       MMVII_INTERNAL_ASSERT_User(mPhProj.DPMulTieP().DirOutIsInit(),eTyUEr::eUnClassedError,"MulTieP In w/o Out");
    }
    tNameSet aSetIn;
    tNameSet aSetOut;
    for (const auto & [aNameIn,aNameOut] :  aDico)
    {
       bool  hasTieP =    isInitMTP
                       && ExistFile(mPhProj.DPMulTieP().FullDirIn()+ mPhProj.NameMultipleTieP(aNameIn));
       int aNbTieP = isInitMTP ? -1 : 0;
       if (hasTieP)
       {
           cVecTiePMul  aVTPM("toto");
           mPhProj.ReadMultipleTieP(aVTPM,aNameIn);
           aVTPM.mNameIm = aNameOut;
           aNbTieP = aVTPM.mVecTPM.size();
           mPhProj.SaveMultipleTieP(aVTPM,aNameOut);
       }
       else
       {
       }

       if (aNbTieP >= mNbMinTieP)
       {
           aSetIn.Add(aNameIn);
           aSetOut.Add(aNameOut);
       }
    }
    SaveInFile(aSetIn ,mNameFilesListIm.at(0));
    SaveInFile(aSetOut,mNameFilesListIm.at(1));

    if (mPhProj.DPGndPt2D().DirInIsInit())
    {
       MMVII_INTERNAL_ASSERT_User(mPhProj.DPGndPt2D().DirOutIsInit(),eTyUEr::eUnClassedError,"Measure In w/o Out");
       for (const auto & aPair :  aDico)
       {
           if (ExistFile(mPhProj.NameMeasureGCPIm(aPair.first,true)))
           {
              cSetMesPtOf1Im  aSMes = mPhProj.LoadMeasureIm(aPair.first);
              aSMes.SetNameIm(aPair.second);
              mPhProj.SaveMeasureIm(aSMes);
           }
       }
       mPhProj.CpGCP();
    }

    return EXIT_SUCCESS;
}

tMMVII_UnikPApli Alloc_DicoRename(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppli_DicoRename(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpecDicoRename
(
    "UtiDicoRename",
    Alloc_DicoRename,
    "This command create a dictionnary after parsing a file, can be used for renaming",
    {eApF::Project},
    {eApDT::FileSys},
    {eApDT::FileSys},
    __FILE__
);



/* ==================================================== */
/*                                                      */
/*          cAppli_Rename                               */
/*                                                      */
/* ==================================================== */


class cAppli_Rename : public cMMVII_Appli
{
     public :
        cAppli_Rename(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli &);  ///< constructor
        int Exe() override;                                             ///< execute action

        int OldExe();
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override; ///< return spec of  mandatory args
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override; ///< return spec of optional args
     protected :
     private :

        std::string ComputeReplace(const std::string & aNameIn) const;

        std::vector<std::string>  Samples() const override;

        void TestSet(const std::string & aName);

        std::string               mPatternGlob; // in simple case, same pat for sel and replace
        std::string               mPatternRepl; // with sub-dir we may need diff pat

        std::string               mSubst;
        std::vector<std::string>  mArithmReplace;
        bool                      mDoReplace;

        std::set<std::string>     mSetOut;
        bool                     mByLink;  /// if true, instead of rename, create a link
        bool                      mShow;  ///< Show msg of replace

};



cCollecSpecArg2007 & cAppli_Rename::ArgObl(cCollecSpecArg2007 & anArgObl)
{
   return anArgObl
            << Arg2007(mPatternGlob,"Pattern of file to replace",{{eTA2007::MPatFile,"0"},{eTA2007::FileDirProj}})
            << Arg2007(mSubst,"Pattern of substituion")
;
}

cCollecSpecArg2007 & cAppli_Rename::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
   return anArgOpt
            << AOpt2007(mDoReplace,"DoReplace","do the replacement ",{{eTA2007::HDV}})
            << AOpt2007(mPatternRepl,"PatRepl","Pattern 4 replace, when != Pattern glob")
            << AOpt2007(mArithmReplace,"AR","arthim repacement like [+,33,2,4] to add 33 to second expr and put on 4 digt ",{{eTA2007::ISizeV,"[3,4]"}})
            << AOpt2007(mByLink,"ByLink","If true create a link instead of moving",{eTA2007::HDV})
            << AOpt2007(mShow,"Show","Show detail of replacment, default = ! DoReplace")
    ;
}


cAppli_Rename::cAppli_Rename(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
  cMMVII_Appli (aVArgs,aSpec),
  mDoReplace   (false),
  mByLink      (false)

{
}

void cAppli_Rename::TestSet(const std::string & aNameOut)
{
    if (BoolFind(mSetOut,aNameOut))
    {
        MMVII_UnclasseUsEr("Proposed replacement would create a conflict");
    }
    mSetOut.insert(aNameOut);
}

std::vector<std::string>  cAppli_Rename::Samples() const
{
  return {
             "MMVII UtiRename \"948_(.*).JPG\" \"\\$&\" AR=[-,1,1] DoReplace=true",
             "MMVII UtiRename  'CCAM_(.*)' 'CCAMBefore_$1'  FFI0=[0000,0030]  PatFFI0=['CCAM_.*_(.*).tif','$1'] DoReplace=true"
         };
}




std::string cAppli_Rename::ComputeReplace(const std::string & aStrIn0) const
{
    //StdOut() << "aStrIn0aStrIn0=[" << aStrIn0 << "]\n";

   std::regex aPat(mPatternRepl);

   std::string aStrIn = aStrIn0;
   if (IsInit(&mArithmReplace))
      aStrIn =  ArithmReplace(aStrIn0, aPat,mArithmReplace);

  //StdOut() << "P="<< mPatternRepl  << " Sub=" << mSubst  << " Str=" << aStrIn << "\n";
   aStrIn =  ReplacePattern(mPatternRepl,mSubst,aStrIn);

   return aStrIn;
}


int cAppli_Rename::Exe()
{
    SetIfNotInit(mShow,!mDoReplace);

    // Not the defaut value is not Pattern glob but mPatOfMS, because with single file
    // the pattern has been replaced by the file (used for parallelization)
    if (!IsInit(&mPatternRepl))
       mPatternRepl = FileOfPath(mPatOfMS[0],false);

    // Compute the map Replace <= [Init1,Init2...]
    bool gotAmbig = false;
    std::map<std::string,std::vector<std::string>> aMapTransfo;
    for (const auto & aStrIn0 : VectMainSet(0))
    {
        std::string aReplace = ComputeReplace(aStrIn0);
        aMapTransfo[aReplace].push_back(aStrIn0);
        if ( aMapTransfo[aReplace].size()>1)
            gotAmbig = true;
        if (mShow)
           StdOut()  << "STR IN=" << aStrIn0 << " => " << aReplace << "\n";
    }

    // is there was any replacement coming from multiple input => error
    if (gotAmbig)
    {
        StdOut() << "=== BAD REPLACE FOR ===============\n";
        for (const auto & [aOut,aIn]: aMapTransfo)
        {
            if (aOut.size()>1)
            {
                StdOut() << " * "<< aOut << " <=== " << aIn << "\n";
            }
        }
        MMVII_UnclasseUsEr("Renaming woul lead to lost file");
    }

    if (mDoReplace)
    {
        for (const auto & [aOut,aVecIn]: aMapTransfo)
        {
            std::string aIn= DirProject() + aVecIn.at(0);
            std::string aDirOut = DirOfPath(aOut,false);
            if ( ExistFile(aDirOut))
                CreateDirectories(aDirOut,false);
          //  StdOut() << "DIOOOO=" << aDirOut  << " " << ExistFile(aDirOut) << "\n";
            if (mByLink)
            {
                CreateLink(aIn,aOut);
            }
            else
            {
                 RenameFiles(aIn,aOut);
            }

        }
    }

    return EXIT_SUCCESS;
}


#if(0)
int cAppli_Rename::OldExe()
{

    std::set<std::string> aSetStr;
    //  StdOut() <<  "============= Proposed replacement  ====== " << std::endl;

    std::vector<std::pair<std::string,std::string>  > aVInOut;

    SetIfNotInit(mShow,!mDoReplace);
    //StdOut() << "RRR " << __LINE__ << "\n";

    std::string aDirLink;
    if (! IsInit(&mPatternRepl))
    {   if (mByLink)
        {
            SplitDirAndFile(aDirLink,mPatternRepl,mPatternGlob,false);

//           mPatternRepl = FileOfPath(mPatternGlob,false);
        }
        else
        {
           mPatternRepl = mPatternGlob;
        }
    }
    std::regex aPat(mPatternRepl);

    //StdOut() << "RRR " << __LINE__ << "\n";

    for (const auto & aStrIn0 : VectMainSet(0))
    {
         //StdOut() << "aStrIn0aStrIn0=[" << aStrIn0 << "]\n";
        std::string aStrIn = aStrIn0;
        if (IsInit(&mArithmReplace))
        {
             std::string anOp = mArithmReplace[0];
             int aOffset = cStrIO<int>::FromStr(mArithmReplace[1]);
             int aKExpr = cStrIO<int>::FromStr(mArithmReplace[2]);

             std::smatch aBoundMatch;
             bool aGotMatch = std::regex_search(aStrIn0, aBoundMatch, aPat);
             Fake4ReleaseUseIt(aGotMatch);
             MMVII_INTERNAL_ASSERT_tiny(aGotMatch,"cCRegex::BoundsMatch no match");

             if ((aKExpr<0)||(aKExpr >= (int) aBoundMatch.size()))
             {
                  MMVII_UnclasseUsEr("Num of expr incompatible with pattern : " + mArithmReplace[2]);
             }

             // auto aMatch  = aBoundMatch[aKExpr];

             std::string aStrNumIn = aBoundMatch[aKExpr];
             int aNum    = cStrIO<int>::FromStr(aStrNumIn);

             if (anOp=="+")
                aNum += aOffset;
             else if (anOp=="-")
                aNum -= aOffset;
             else if (anOp=="%")
                aNum %= aOffset;
             else
             {
                MMVII_UnclasseUsEr("Bad operand in arithmetic : " + mArithmReplace[0]);
             }

             int aNbDig = (mArithmReplace.size()> 3) ? cStrIO<int>::FromStr(mArithmReplace[3]) : aStrNumIn.size();
             std::string aStrNumOut = ToStr(aNum,aNbDig);

             aStrIn.replace(aBoundMatch.position(aKExpr),aBoundMatch.length(aKExpr),aStrNumOut);
        }
        if (0)
        {
            StdOut() << "P=" << mPatternRepl << " S=" << mSubst << " I=" << aStrIn << "\n";
        }
        std::string aStrOut =  ReplacePattern(mPatternRepl,mSubst,aStrIn);
        // StdOut() << "[" << aStrIn0  << "] ";
        if (IsInit(&mArithmReplace))
        {
            if (mShow)
               StdOut() << " AR==> [" << aStrIn  << "] ";
        }

        if (mShow)
            StdOut() << " ==> [" << aStrOut  << "]  " << std::endl;

        // TestSet(aStrIn0);
        TestSet(aStrOut);
        aVInOut.push_back(std::pair<std::string,std::string>(aStrIn0,aStrOut));
    }

    for (const auto & aPair : aVInOut)
    {
       // auto [aStrIn0,aStrOut] = aPair;
       auto aStrOut = aPair.second;
       if (ExistFile(aStrOut) && (! BoolFind(mSetOut,aStrOut)))
       {
           MMVII_UnclasseUsEr("File already exist");
       }
    }
    //StdOut() << " NbFiles= " << aVInOut.size() << "\n";

    std::string aPrefTmp = "MMVII_Tmp_Replace_"+ PrefixGMA() + "_";

    if (mByLink)
    {
       for (const auto &  [aStrIn0,aStrOut]  : aVInOut)
       {
        //   StdOut()  << " LLLnk " << aDirLink+aStrIn0 << " " << aStrOut << "\n";
           if (mDoReplace)
               CreateLink(aDirLink+aStrIn0,aStrOut);
       }
    }
    else
    {
        if (mDoReplace)
        {
            // In case "input" intersect "outout", put first "input" in "tmp" file,
            for (const auto & aPair : aVInOut)
            {
                auto [aStrIn0,aStrOut] = aPair;
                StdOut() << "mv " << aStrIn0  << " " << aPrefTmp+aStrIn0  << std::endl;
                RenameFiles(aStrIn0,aPrefTmp+aStrIn0);
            }
            // the put, safely, "tmp" in "output"
            for (const auto & aPair : aVInOut)
            {
                auto [aStrIn0,aStrOut] = aPair;
                StdOut() << "mv " << aPrefTmp+ aStrIn0  << " " << aStrOut  << std::endl;
                RenameFiles(aPrefTmp+aStrIn0,aStrOut);
            }
        }
    }

    return EXIT_SUCCESS;
}
#endif


tMMVII_UnikPApli Alloc_Rename(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppli_Rename(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpecRename
(
    "UtiRename",
    Alloc_Rename,
    "This command renames files using regexpr and eventually arithmetic",
    {eApF::Project},
    //  {eApF::ManMMVII, eApF::Project},  JOE ?  j'ai enleve eApF::ManMMVI, je sais plus qui l'a mis
    {eApDT::FileSys},
    {eApDT::FileSys},
    __FILE__
);

}
