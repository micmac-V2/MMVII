#include "cMMVII_Appli.h"
#include "MMVII_Tpl_ElemStrToVal.h"
#include "MMVII_Sys.h"
#include "MMVII_DeclareCste.h"
#include "MMVII_2Include_Serial_Tpl.h"
#include "../Serial/Serial.h"


#include <functional>
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


Color Color::command;
Color Color::argument;
Color Color::title;
Color Color::success;
Color Color::error;
Color Color::warning;
Color Color::descr;
Color Color::end;

bool Color::mNoColor = true;

static void SetColors(eModeHelpColor aMode)
{
    static const Color Reset       ("\033[0m");

    static const Color DarkRed     ("\033[31m");
    static const Color DarkGreen   ("\033[32m");
    static const Color DarkYellow  ("\033[33m");
    static const Color DarkBlue    ("\033[34m");
    static const Color DarkMagenta ("\033[35m");
    static const Color DarkCyan    ("\033[36m");
    static const Color LightGray   ("\033[37m");

    static const Color DarkGray    ("\033[90m");
    static const Color Red         ("\033[91m");
    static const Color Green       ("\033[92m");
    static const Color Yellow      ("\033[93m");
    static const Color Blue        ("\033[94m");
    static const Color Magenta     ("\033[95m");
    static const Color Cyan        ("\033[96m");
    static const Color White       ("\033[97m");

    if (aMode == eModeHelpColor::None)
    {
        Color::Off();
        return;
    }
    Color::On();
    if (aMode == eModeHelpColor::Light)
    {
        Color::command = DarkBlue;
        Color::argument = DarkMagenta;
        Color::descr = DarkGreen;
        Color::title = DarkCyan;
        Color::success = DarkGreen;
        Color::error = DarkRed;
        Color::warning = DarkMagenta;
        Color::end = Reset;
    }
    else
    {
        Color::command = DarkGreen;
        Color::argument = DarkGreen;
        Color::descr = DarkYellow;
        Color::success = Green;
        Color::title = Blue;
        Color::error = Red;
        Color::warning = DarkYellow;
        Color::end = Reset;
    }
}


/* ==================================================== */
/*                                                      */
/*          Helpers for profile managment               */
/*                                                      */
/* ==================================================== */

extern std::string MMVII_UserConfigDir(bool aSVP);

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

/// Values used when no profile at all can be read, they must be valid for any command
static cParamProfile BuiltInProfile()
{
    cParamProfile aRes;
    aRes.Set("VectSerialMode",eTypeSerial::ecsv);
    aRes.Set("TaggedSerialMode",eTypeSerial::exml);
    aRes.Set("HelpColorMode",eModeHelpColor::Dark);
    return aRes;
}

static cParamProfile SyncUserProfile(const std::string & aUserFile,const std::string & aDefaultFile)
{
    cParamProfile aDefaultProfile;
    if (! ExistFile(aDefaultFile))
    {                 // Shouldn't happen ...
        aDefaultProfile = BuiltInProfile();
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
    mDirUserProfile = MMVII_UserConfigDir(SVP::Yes);
    if (mDirUserProfile.empty())   // no user config dir, see InitProfile
       return "Default";
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
    mDirUserProfile = MMVII_UserConfigDir(SVP::Yes);
    if (mDirUserProfile.empty())
    {
        //  The user configuration directory cannot even be named, typically because HOME is
        //  not defined.  Nothing can be read from or written to a profile, but the command
        //  itself must still run, with the built-in values.
        MMVII_USER_WARNING("Cannot locate the user configuration directory, running with default profile values");
        mUserProfile = "";
        mParamProfile = BuiltInProfile();
    }
    else
    {
        MakeNameDir(mDirUserProfile);

        if (! IsInit(&mProfileName)) {
            CreateDirectories(mDirUserProfile,SVP::No);

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
    }
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
    //  When there is no user configuration directory, no profile file can be enumerated,
    //  only the default one exists, in memory
    if (! aDirUserProfile.empty())
    {
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


class cValidator
{
public:
    template <typename T> static cValidator ForType()
    {
        std::string aAllowed;
        if constexpr (std::is_enum_v<T>) {
            aAllowed = StrAllVall<T>();
        }
        return cValidator
        (
            [](const std::string & aVal)
            {
                try
                {
                    cStrIO<T>::FromStr(aVal,true);
                    return true;
                }
                catch (const StrIOException &)
                {
                    return false;
                }
            },
            cStrIO<T>::msNameType(),
            aAllowed
        );
    }

    bool operator()(const std::string & aVal) const
    {
        return mValidate(aVal);
    }

    const std::string & TypeName() const { return mTypeName; }
    const std::string & Allowed() const { return mAllowedValues; }

private:
    cValidator
    (
        const std::function<bool(const std::string &)> & aValidate,
        const std::string & aTypeName,
        const std::string & aAllowedValues
    ) :
        mTypeName (aTypeName),
        mAllowedValues (aAllowedValues),
        mValidate (aValidate)
    {
    }

    std::string mTypeName;
    std::string mAllowedValues;
    std::function<bool(const std::string &)> mValidate;
};

using tMapProfileValidator = std::map<std::string,cValidator>;

static const tMapProfileValidator & MapProfileKeyValidator()
{
    static const tMapProfileValidator TheMap
    {
        {"HelpColorMode",    cValidator::ForType<eModeHelpColor>()},
        {"NbProcMax",        cValidator::ForType<int>()},
        {"PdfOpen",          cValidator::ForType<std::string>()},
        {"TaggedSerialMode", cValidator::ForType<eTypeSerial>()},
        {"VectSerialMode",   cValidator::ForType<eTypeSerial>()}
    };
    return TheMap;
}

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
        << AOpt2007(mKeyVal,"KeyVal","Set Value to Key: [Key,Value]",{{eTA2007::ISizeV,"[2,2]"}})
        << AOpt2007(mDelKey,"DelKey","Delete Key",{eTA2007::ProfileKey})
   ;
}



int cAppli_EditProfile::Exe()
{
    //  Startup can go on with default values when there is no user configuration directory,
    //  but this command was explicitly called to show or modify a profile : here it is an error
    MMVII_INTERNAL_ASSERT_User
    (
        ! mDirUserProfile.empty(),
        eTyUEr::eUnClassedError,
        "Cannot locate the user configuration directory, check the environment (HOME ...)"
    );

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

        for (const auto & [aKey,aValidator] : MapProfileKeyValidator())
        {
            (void) aValidator;
            std::string aVal ="";
            if ( mParamProfile.HasKey(aKey))
                aVal =  mParamProfile.Get(aKey,std::string("<Missing>"));
            StdOut() << "  " << aKey << "= \"" <<  aVal  << "\"\n";
        }
        return EXIT_SUCCESS;
    }

    const std::string aNameProfile = aCurrent ? mCurrent : ProfileName();
    CheckProfileName(aNameProfile);
    const std::string aFileParam = mDirUserProfile + TheUserProfilePrefix + aNameProfile + ".xml";
    cParamProfile aParam = SyncUserProfile(aFileParam,mDirLocalParameters+TheDefaultProfileName);

    if (IsInit(&mKeyVal))
    {
        const auto& aKey = mKeyVal[0];
        const auto & aMapValidator = MapProfileKeyValidator();
        const auto aItValidator = aMapValidator.find(aKey);

        if (aItValidator != aMapValidator.end())
        {
            const auto& aVal = mKeyVal[1];

            auto& aValidator = aItValidator->second;
            MMVII_INTERNAL_ASSERT_User
            (
                aItValidator->second(aVal),
                eTyUEr::eBadOptParam,
                "Invalid value [" + aVal + "] for profile key [" + aKey
                    + "], expecting " + ( aValidator.Allowed().empty() ? ("type " + aValidator.TypeName()) : aValidator.Allowed() )
            );
            StdOut() << "Setting profile key '" << aKey << "' to value '" << aVal << "'\n";
            aParam.Set(aKey, aVal);
        }
        else
        {
            MMVII_USER_WARNING("Unknown profile key: " + aKey);
            StdOut() << "--------- Allowed values -------------\n";
            for (const auto & [aKeyAllowed,aValidator] : aMapValidator)
            {
                (void) aValidator;
                StdOut() << "  * " << aKeyAllowed << "\n";
            }
        }
    }
    if (IsInit(&mDelKey))
    {
        if (aParam.HasKey(mDelKey))
        {
           StdOut() << "Deleting key '" << mDelKey << "'\n";
           aParam.Remove(mDelKey);
        }
        else
        {
            MMVII_USER_WARNING("Key to remove doesn't exist : " + mDelKey);
        }
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

} // namespace MMVII

