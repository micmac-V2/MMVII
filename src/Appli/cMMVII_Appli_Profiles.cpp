#include "cMMVII_Appli.h"
#include "MMVII_Sys.h"
#include "MMVII_DeclareCste.h"
#include "MMVII_2Include_Serial_Tpl.h"


#include <thread>
#include <mutex>

namespace MMVII
{


struct cSpecifProfileUserMMVII
{
     public :
         std::string mNameProfile;
};

void AddData(const cAuxAr2007 & anAux,cSpecifProfileUserMMVII & aSpec)
{
     AddData(cAuxAr2007("NameProfile",anAux),aSpec.mNameProfile);
}

void AddData(const cAuxAr2007 & anAux,cParamProfile & aProfile)
{
     AddData(cAuxAr2007("UserName",anAux),aProfile.mUserName);
     AddData(cAuxAr2007("NbProcMax",anAux),aProfile.mNbProcMax);
     EnumAddData(anAux,aProfile.mTaggedDefSerial,"TaggedSerialMode");
     EnumAddData(anAux,aProfile.mVectDefSerial,"VectSerialMode");
}


/* ==================================================== */
/*                                                      */
/*          cMMVII_Appli                                */
/*                                                      */
/* ==================================================== */

cParamProfile cMMVII_Appli::mParamProfile;
std::string cMMVII_Appli::mDirProfileUsage;
std::string cMMVII_Appli::mProfileUsage;


std::string TheNameFileCurentProfile =  "MMVII-CurentProfile.xml";
std::string NameFileUseOfProfile =  "MMVII-UserOfProfile.xml";
std::string DefaultNameFileCurentProfile()
{
    return  "Default-" + TheNameFileCurentProfile;
}


cParamProfile::cParamProfile() :
   mUserName        ("Unknown"),
   mNbProcMax       (1000),
   mTaggedDefSerial (eTypeSerial::exml),
   mVectDefSerial   (eTypeSerial::ecsv)

{
}

void cMMVII_Appli::InitProfile()
{

  // ========================================================================
  // ========================  HANDLING PROFILE USER ETC ... ================
  // ========================================================================


  //  part of code that was used to initialize "at hand", soon will be obsolete...
  if (0)
  {
      StdOut() << "NO USEERRRRRRRRRRRRRR " << std::endl; getchar();

      mParamProfile.mUserName = "Unknown";
      mParamProfile.mNbProcMax = 1000;
      mParamProfile.mVectDefSerial = eTypeSerial::ecsv;
      mParamProfile.mTaggedDefSerial = eTypeSerial::ejson;
      mVectNameDefSerial = E2Str(mParamProfile.mVectDefSerial);
      mTaggedNameDefSerial = E2Str(mParamProfile.mTaggedDefSerial);
      return;
  }

  /*  Compute the name of file containing the profile of user;  this profile is
   *
   *     - "MMVII-CurentPofile.xml" if this file exists, to allow tuning by user
   *     - "Default-MMVII-CurentPofile.xml" if it does not exist, this is the file shared on github
   */

  // if the default file  does not exist, we are probably the first time, or in reinit step because
  // directory has been purged, we create a file containing  the default profile
  if (! ExistFile(mDirLocalParameters+DefaultNameFileCurentProfile()))
  {
        cSpecifProfileUserMMVII  aSpec;
        aSpec.mNameProfile = "Default";
        SaveInFile(aSpec,mDirLocalParameters+DefaultNameFileCurentProfile());
  }

  // we set NameFileCurentProfile  to its default or not,
  // StdOut() <<  "EXXXX " << ExistFile(mDirLocalParameters+TheNameFileCurentProfile)
  //         << " " << mDirLocalParameters+TheNameFileCurentProfile << "\n";
  if (! ExistFile(mDirLocalParameters+TheNameFileCurentProfile))
  {
      TheNameFileCurentProfile = DefaultNameFileCurentProfile();
  }

  /**  Compute the "usage" store in the profile,
   *   init the variable  "mProfileUsage"  and  "mDirProfileUsage"
   */
  {
      cSpecifProfileUserMMVII aSpec;

      ReadFromFile(aSpec,mDirLocalParameters+TheNameFileCurentProfile);
      mProfileUsage = aSpec.mNameProfile;
      mDirProfileUsage =  mDirLocalParameters + mProfileUsage + StringDirSeparator();
  }

 // StdOut() << "IIiiiiiiiiiiiiiiiiiINiProfileEEEEEEEEE\n";
  /**  if file containing users profile does not exist, we create some default one */
  if (! ExistFile(mDirProfileUsage+NameFileUseOfProfile))
  {
      CreateDirectories(mDirProfileUsage,false);

      mParamProfile.mUserName = "Unknown";
      mParamProfile.mNbProcMax = 1000;
      mParamProfile.mVectDefSerial = eTypeSerial::ecsv;
      mParamProfile.mTaggedDefSerial = eTypeSerial::exml;
      SaveInFile(mParamProfile,mDirProfileUsage+NameFileUseOfProfile);
  }
  ReadFromFile(mParamProfile,mDirProfileUsage+NameFileUseOfProfile);
  mVectNameDefSerial = E2Str(mParamProfile.mVectDefSerial);
  mTaggedNameDefSerial = E2Str(mParamProfile.mTaggedDefSerial);

}

const  std::string & cMMVII_Appli::UserName() {return mParamProfile.mUserName;}
const  std::string & cMMVII_Appli::DirProfileUsage() {return mDirProfileUsage;}

eTypeSerial cMMVII_Appli::VectDefSerial() const
{
    CurrentAppli(); // as member is static assure init was done
    return mParamProfile.mVectDefSerial;
}
eTypeSerial cMMVII_Appli::TaggedDefSerial() const
{
    CurrentAppli(); // as member is static assure init was done
    return mParamProfile.mTaggedDefSerial;
}

const std::string & cMMVII_Appli::VectNameDefSerial   () const
{
    CurrentAppli(); // as member is static assure init was done
    return mVectNameDefSerial;
}
const std::string & cMMVII_Appli::TaggedNameDefSerial   () const
{
    CurrentAppli(); // as member is static assure init was done
    return mTaggedNameDefSerial;
}

const std::string & GlobVectNameDefSerial() {return cMMVII_Appli::CurrentAppli().VectNameDefSerial();}
const std::string & GlobTaggedNameDefSerial() {return cMMVII_Appli::CurrentAppli().TaggedNameDefSerial();}



/* ==================================================== */
/*                                                      */
/*          cAppli_EditProfile                          */
/*                                                      */
/* ==================================================== */


class cAppli_EditProfile : public cMMVII_Appli
{
     public :
        cAppli_EditProfile(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli &);
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;

     private :

         bool            mSetUser;
         cParamProfile   mModifParam;
};

cAppli_EditProfile::cAppli_EditProfile(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
    cMMVII_Appli  (aVArgs,aSpec),
    mSetUser      (false)
{

}

cCollecSpecArg2007 & cAppli_EditProfile::ArgObl(cCollecSpecArg2007 & anArgObl)
{
   return
      anArgObl
         << Arg2007(mModifParam.mUserName,"Name of user")
      ;
}

cCollecSpecArg2007 & cAppli_EditProfile::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
   return
      anArgOpt
        << AOpt2007(mSetUser,"SetUser","Set this user as the current user of MMVII")
        << AOpt2007(mModifParam.mNbProcMax,"DefNbProc","Default number of processor")
        << AOpt2007(mModifParam.mTaggedDefSerial,"ModeSerial","Mode for serialization",{AC_ListVal<eTypeSerial>()})
    ;
}

int cAppli_EditProfile::Exe()
{
    cParamProfile aParam;

    std::string aDirParam = mDirLocalParameters + mModifParam.mUserName ;
    std::string aFileParam = aDirParam + StringDirSeparator()+ NameFileUseOfProfile;

    if (ExistFile(aFileParam))
    {
        ReadFromFile(aParam,aFileParam);
    }
    else
    {
        CreateDirectories(aDirParam,false);
    }

    aParam.mUserName = mModifParam.mUserName;
    if (IsInit(&mModifParam.mNbProcMax))
        aParam.mNbProcMax = mModifParam.mNbProcMax;
    if (IsInit(&mModifParam.mTaggedDefSerial))
        aParam.mTaggedDefSerial = mModifParam.mTaggedDefSerial;

    SaveInFile(aParam,aFileParam);

    if (mSetUser)
    {
        std::string aFileSpec = mDirLocalParameters + TheNameFileCurentProfile;
       // StdOut() << "FSSS=" << aFileSpec << " UN=" << mUserName << "\n";
        cSpecifProfileUserMMVII aSpec;
        ReadFromFileWithDef(aSpec,aFileSpec);

        aSpec.mNameProfile = mModifParam.mUserName;

        SaveInFile(aSpec,aFileSpec);
    }
    else
    {
        StdOut() << " Curent user="  << UserName()   << ",not modified\n";
    }

    return EXIT_SUCCESS;
}

tMMVII_UnikPApli Alloc_EditProfile(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppli_EditProfile(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpecEditProfile
(
     "EditProfile",
      Alloc_EditProfile,
      "This command is used to edit a user's profile ",
      {eApF::Project},
      {eApDT::Console,eApDT::Xml},
      {eApDT::Xml},
      __FILE__
);



};

