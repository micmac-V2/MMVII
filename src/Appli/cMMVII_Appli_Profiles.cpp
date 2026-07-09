#include "cMMVII_Appli.h"
#include "MMVII_Tpl_ElemStrToVal.h"
#include "MMVII_Sys.h"
#include "MMVII_DeclareCste.h"
#include "MMVII_2Include_Serial_Tpl.h"
#include "../Serial/Serial.h"


#include <thread>
#include <mutex>

namespace MMVII
{

/* ==================================================== */
/*                                                      */
/*          Colors for Help                             */
/*                                                      */
/* ==================================================== */

/// Color set used by command help and command list output
enum class eModeHelpColor
{
    Dark,      ///< Color set for dark terminal background
    Light,     ///< Color set for light terminal background
    None,      ///< Disable terminal colorization
    eNbVals    ///< Tag for number of value
};

template<> cE2Str<eModeHelpColor>::tMapE2Str cE2Str<eModeHelpColor>::mE2S
    {
        {eModeHelpColor::Dark,"Dark"},
        {eModeHelpColor::Light,"Light"},
        {eModeHelpColor::None,"None"}
    };

MACRO_INSTANTIATE_STRIO_ENUM(eModeHelpColor,"ModeHelpColor")


// TODOCM: Supprimer std:: de std::vector dans HElp ?
// TODOCM: SetColors() ou executer ? Apres Profile ?  Attention au InitProfile dans GenerateHelp Et dans main() pour help general ?
// TODOCM: profile : possibilite check values ...

Col Col::Reset       ("\033[0m");

Col Col::DarkRed     ("\033[31m");
Col Col::DarkGreen   ("\033[32m");
Col Col::DarkYellow  ("\033[33m");
Col Col::DarkBlue    ("\033[34m");
Col Col::DarkMagenta ("\033[35m");
Col Col::DarkCyan    ("\033[36m");
Col Col::LightGray   ("\033[37m");

Col Col::DarkGray    ("\033[90m");
Col Col::Red         ("\033[91m");
Col Col::Green       ("\033[92m");
Col Col::Yellow      ("\033[93m");
Col Col::Blue        ("\033[94m");
Col Col::Magenta     ("\033[95m");
Col Col::Cyan        ("\033[96m");
Col Col::White       ("\033[97m");

Col Col::command;
Col Col::argument;
Col Col::title;
Col Col::error;
Col Col::warning;
Col Col::descr;
Col Col::end;

static void SetColors(eModeHelpColor aMode)
{
    if (aMode == eModeHelpColor::None)
    {
        Col::command = Col();
        Col::argument =Col();
        Col::descr = Col();
        Col::title = Col();
        Col::error = Col();
        Col::warning = Col();
        Col::end = Col();
    }
    else if (aMode == eModeHelpColor::Light)
    {
        Col::command = Col::DarkBlue;
        Col::argument = Col::DarkMagenta;
        Col::descr = Col::DarkGreen;
        Col::title = Col::DarkCyan;
        Col::error = Col::DarkRed;
        Col::warning = Col::DarkMagenta;
        Col::end = Col::Reset;
    }
    else
    {
        Col::command = Col::DarkGreen;
        Col::argument = Col::DarkGreen;
        Col::descr = Col::DarkYellow;
        Col::title = Col::Blue;
        Col::error = Col::Red;
        Col::warning = Col::DarkYellow;
        Col::end = Col::Reset;
    }
}


/* ==================================================== */
/*                                                      */
/*          Helpers for profile managment               */
/*                                                      */
/* ==================================================== */

extern std::string MMVII_UserConfigDir();

struct cSelectedProfile
{
    std::string mNameProfile;
};

void AddData(const cAuxAr2007 & anAux,cSelectedProfile & aProfile)
{
    AddData(cAuxAr2007("NameProfile",anAux),aProfile.mNameProfile);
}

static void CheckProfileName(const std::string & aName)
{
    MMVII_INTERNAL_ASSERT_User
    (
        (!aName.empty()) && (aName==ToStandardStringIdent(aName)),
        eTyUEr::eBadOptParam,
        "Profile name must contain only letters, digits, underscore or hyphen"
    );
}

static cParamProfile SyncUserProfile(const std::string & aUserFile,const std::string & aDefaultFile)
{
    cParamProfile aDefaultProfile;
    if (! ExistFile(aDefaultFile))
    {                 // Shouldn't happen ...
        aDefaultProfile.Set("VectSerialMode",eTypeSerial::ecsv);
        aDefaultProfile.Set("TaggedSerialMode",eTypeSerial::exml);
        aDefaultProfile.Set("HelpColorMode",eModeHelpColor::Dark);
        SaveInFile(aDefaultProfile,aDefaultFile);          //  ... create one just in case
    }
    ReadFromFile(aDefaultProfile,aDefaultFile);
    bool aModified = false;
    cParamProfile aUserProfile;
    if (ExistFile(aUserFile))
    {
        ReadFromFile(aUserProfile,aUserFile);
        for (const auto& [aKey, aVal] : aDefaultProfile)
        {
            if (! aUserProfile.HasKey(aKey))
            {
                aUserProfile.Set(aKey,aVal);
                aModified = true;
            }
        }
    } else {
        aUserProfile = aDefaultProfile;
        aModified = true;
    }
    if (aModified)
    {
        SaveInFile(aUserProfile,aUserFile);
    }
    return aUserProfile;
}


void AddData(const cAuxAr2007 & anAux,cParamProfile & aProfile)
{
    aProfile.AddData(anAux);
}

cParamProfile::cParamProfile()
{
}


void cParamProfile::AddData(const cAuxAr2007 &anAux)
{
    StdMapAddData(anAux,mMap);
}


/* ==================================================== */
/*                                                      */
/*          cMMVII_Appli                                */
/*                                                      */
/* ==================================================== */

cParamProfile cMMVII_Appli::mParamProfile;
std::string cMMVII_Appli::mDirUserProfile;
std::string cMMVII_Appli::mProfileName;
std::string cMMVII_Appli::mUserProfile;


const std::string TheDefaultProfileName = "MMVII-Default-Profile.xml";
const std::string TheSelectedProfileName = "MMVII-Current-Profile.xml";
const std::string TheUserProfilePrefix = "MMVII-profile-";

std::string cMMVII_Appli::GetProfileName()
{
    mDirUserProfile = MMVII_UserConfigDir();
    MakeNameDir(mDirUserProfile);
//    CreateDirectories(mDirUserProfile,false);
    const std::string aSelectedProfile = mDirUserProfile + TheSelectedProfileName;
    cSelectedProfile aSelection;
    if (!ExistFile(aSelectedProfile))
    {
        aSelection.mNameProfile = "Default";
    } else {
        ReadFromFile(aSelection,aSelectedProfile);
    }
    return aSelection.mNameProfile;
}

void cMMVII_Appli::InitProfile()
{
    mDirUserProfile = MMVII_UserConfigDir();
    MakeNameDir(mDirUserProfile);

    if (! IsInit(&mProfileName)) {
        CreateDirectories(mDirUserProfile,false);

        const std::string aSelectedProfile = mDirUserProfile + TheSelectedProfileName;
        if (!ExistFile(aSelectedProfile))
        {
            cSelectedProfile aSelection;
            aSelection.mNameProfile = "Default";
            SaveInFile(aSelection,aSelectedProfile);
        }
        cSelectedProfile aSelection;
        ReadFromFile(aSelection,aSelectedProfile);
        mProfileName = aSelection.mNameProfile;
    }
    CheckProfileName(mProfileName);

    mUserProfile = mDirUserProfile + TheUserProfilePrefix + mProfileName + ".xml";
    const std::string aDefaultProfile = mDirLocalParameters + TheDefaultProfileName;
    mParamProfile = SyncUserProfile(mUserProfile,aDefaultProfile);
    mVectNameDefSerial = mParamProfile.Get("VectSerialMode",ToS(eTypeSerial::ecsv));
    mTaggedNameDefSerial = mParamProfile.Get("TaggedSerialMode",ToS(eTypeSerial::exml));
    if (mParamProfile.HasKey("HelpColorMode")) {
        SetColors(mParamProfile.Get("HelpColorMode",eModeHelpColor::Dark));
    } else {
        SetColors(eModeHelpColor::Dark);
    }
}

const std::string & cMMVII_Appli::ProfileName() {return mProfileName;}
const std::string & cMMVII_Appli::DirUserProfile() {return mDirUserProfile;}

eTypeSerial cMMVII_Appli::VectDefSerial() const
{
    CurrentAppli(); // as member is static assure init was done
    return mParamProfile.Get("VectSerialMode",eTypeSerial::ecsv);
}
eTypeSerial cMMVII_Appli::TaggedDefSerial() const
{
    CurrentAppli(); // as member is static assure init was done
    return mParamProfile.Get("TaggedSerialMode",eTypeSerial::exml);
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

std::vector<std::string> GlobProfileNames()
{
    auto aDirUserProfile = cMMVII_Appli::DirUserProfile();
    const std::string aSuffix = ".xml";
    std::vector<std::string> aNames;
    for (const auto & aFile : GetFilesFromDir(aDirUserProfile,AllocRegex(TheUserProfilePrefix+".*\\.xml")))
    {
        aNames.push_back
        (
            aFile.substr
            (
                TheUserProfilePrefix.size(),
                aFile.size()-TheUserProfilePrefix.size()-aSuffix.size()
            )
        );
    }
    if (aNames.empty())
    {
        aNames.push_back("Default");
    }
    std::sort(aNames.begin(),aNames.end());
    return aNames;
}


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

         std::string     mCurrent;
         std::vector<std::string>   mKeyVal;
         std::string   mDelKey;
         eModeHelpColor mColor;
};

cAppli_EditProfile::cAppli_EditProfile(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
    cMMVII_Appli  (aVArgs,aSpec)
{

}

cCollecSpecArg2007 & cAppli_EditProfile::ArgObl(cCollecSpecArg2007 & anArgObl)
{
   return anArgObl;
}

cCollecSpecArg2007 & cAppli_EditProfile::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
   return
      anArgOpt
        << AOpt2007(mCurrent,"SetCurrent","Name of the profile to make current",{eTA2007::Profile})
        << AOpt2007(mKeyVal,"KeyVal","Set Value to Key: [Key,Value]")
        << AOpt2007(mDelKey,"DelKey","Delete Key",{eTA2007::ProfileKey})
           << AOpt2007(mColor,"Color","Delete Key")
   ;
}

int cAppli_EditProfile::Exe()
{
    const bool aCurrent = IsInit(&mCurrent);
    const bool aHasModification = aCurrent || IsInit(&mKeyVal) || IsInit(&mDelKey);

    if (!aHasModification)
    {
        const std::vector<std::string> aNames = GlobProfileNames();

        StdOut() << "Profile names:";
        for (const auto & aName : aNames)
            StdOut() << " " << ((aName==ProfileName()) ? "["+aName+"]" : aName);
        StdOut() << "\n";
        StdOut() << "Current profile path is " << mUserProfile << "\n";
        StdOut() << "Profile values:\n";
        for (const auto& [aKey, aVal] : mParamProfile)
            StdOut() << "  " << aKey << ": " << aVal << "\n";
        return EXIT_SUCCESS;
    }

    const std::string aNameProfile = aCurrent ? mCurrent : ProfileName();
    CheckProfileName(aNameProfile);
    const std::string aFileParam = mDirUserProfile + TheUserProfilePrefix + aNameProfile + ".xml";
    cParamProfile aParam = SyncUserProfile(aFileParam,mDirLocalParameters+TheDefaultProfileName);

    if (IsInit(&mKeyVal))
    {
        if (mKeyVal.size() != 2)
            MMVII_UserError(eTyUEr::eBadSize4Vect, "Each KeyVal must have exactly two elements: key and value");
        const auto& aKey = mKeyVal[0];
        const auto& aVal = mKeyVal[1];
        StdOut() << "Setting profile key '" << aKey << "' to value '" << aVal << "'\n";
        aParam.Set(aKey, aVal);
    }
    if (IsInit(&mDelKey))
    {
        StdOut() << "Deleting key '" << mDelKey << "'\n";
        aParam.Remove(mDelKey);
    }
    SaveInFile(aParam,aFileParam);

    if (aCurrent)
    {
        cSelectedProfile aSelection;
        aSelection.mNameProfile = aNameProfile;
        SaveInFile(aSelection,mDirUserProfile+TheSelectedProfileName);
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
      "Edit the current profile ",
      {eApF::Project},
      {eApDT::Console,eApDT::Xml},
      {eApDT::Xml},
      __FILE__
);



};

