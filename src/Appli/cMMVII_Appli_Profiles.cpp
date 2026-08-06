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
        Color::error = DarkRed; // Mpd to have descr != warn
        Color::warning = Red;
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

/**  Thrown by cProfileErrorCatcher, caught in this file, never crosses its interface */
struct cProfileError
{
    std::string mType;
    std::string mMes;
};

static void ProfileErrorHandler(const std::string & aType,const std::string & aMes,const char *,int)
{
    throw cProfileError{aType,aMes};
}

/**  Could the profile not be *reached*, as opposed to not be understood ?  Being unable to
     create the directory or to open the file is an ordinary situation -- a read only home, a
     first run on a system where nothing may be written -- which must stay silent. An
     unusable content is different: the file is there, the user believes their settings
     apply, and they silently do not; that one deserves a warning.
     Built with E2Str so that renaming an enumerator cannot silently break the test.
*/
static bool IsUnreachableProfile(const std::string & aTypeOfError)
{
    return    (aTypeOfError == "UserEr:" + E2Str(eTyUEr::eOpenFile))
           || (aTypeOfError == "UserEr:" + E2Str(eTyUEr::eCreateDir));
}

/**  While such an object is alive, a MMVII error throws a cProfileError instead of aborting.
     Used to make the profile initialization non fatal : whatever fails (unwriteable directory,
     corrupted or truncated file, disk full ...) the command must go on with its default values.
     Reading errors are safe to unwind : the tree of a corrupted file is built in the constructor
     of the archive, before any cAuxAr2007 exists, and cIMakeTreeAr sets its mError flag before
     raising, so the destructors running during the unwinding do not raise in their turn.
*/
class cProfileErrorCatcher
{
   public :
      cProfileErrorCatcher() :
          mPrevHandler (MMVVI_Error)   // may be not the default one, i.e. the python API
      {
          MMVII_SetErrorHandler(ProfileErrorHandler);
      }
      ~cProfileErrorCatcher() {MMVII_SetErrorHandler(mPrevHandler);}
   private :
      cProfileErrorCatcher(const cProfileErrorCatcher &) = delete;
      PtrMMVII_Error_Handler mPrevHandler;
};

/// Why the profile could not be used, empty if it could; kept for the EditProfile message
static std::string TheProfileFailure;

/**  To be called before SaveInFile : the tagged archives write the file in their destructor,
     so a file that cannot be written raises there, where cProfileError cannot be thrown.
     Opening it here, in append mode so that nothing is lost, raises at a catchable place.
*/
static void CheckCanWrite(const std::string & aNameFile)
{
    cMMVII_Ofs aOfs(aNameFile,eFileModeOut::AppendText);
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
        CheckCanWrite(aDefaultFile);
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
        CheckCanWrite(aUserFile);
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
    //  Called for each application of the process, also for those allocated by GenArgsSpec :
    //  once the profile has been diagnosed as unusable, don't restore a directory that would
    //  make the next applications believe it can be used
    if (! TheProfileFailure.empty())
    {
       mDirUserProfile = "";
       return "Default";
    }
    mDirUserProfile = MMVII_UserConfigDir(SVP::Yes);
    if (mDirUserProfile.empty())   // no user config dir, see InitProfile
       return "Default";
    try
    {
        //  Called for every command, even those that never use a profile : a selection file
        //  that cannot be read must not prevent them from running. InitProfile does the diagnostic.
        cProfileErrorCatcher aCatchErrors;

        MakeNameDir(mDirUserProfile);
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
    catch (const cProfileError &)
    {
        return "Default";
    }
}

void cMMVII_Appli::InitProfile()
{
    TheProfileFailure = "";
    //  A profile that cannot even be located or reached is not worth a warning, see below
    bool aProfileIsUnusable = false;
    mDirUserProfile = MMVII_UserConfigDir(SVP::Yes);
    if (mDirUserProfile.empty())
    {
        //  The user configuration directory cannot even be named, typically because HOME is
        //  not defined
        TheProfileFailure = "cannot locate the user configuration directory, check the environment (HOME ...)";
    }
    else
    {
        //  A profile name given on the command line is an explicit user choice : an invalid
        //  one remains a hard error, it is not something to degrade silently
        if (IsInit(&mProfileName))
            CheckProfileName(mProfileName);
        try
        {
            //  Any failure below (unwriteable directory, corrupted file, disk full ...) must
            //  degrade to the built-in values, the command itself has still to run
            cProfileErrorCatcher aCatchErrors;

            MakeNameDir(mDirUserProfile);

            if (! IsInit(&mProfileName)) {
                CreateDirectories(mDirUserProfile);

                const std::string aSelectedProfile = mDirUserProfile + TheSelectedProfileName;
                if (!ExistFile(aSelectedProfile))
                {
                    cSelectedProfile aSelection;
                    aSelection.mNameProfile = "Default";
                    CheckCanWrite(aSelectedProfile);
                    SaveInFile(aSelection,aSelectedProfile);
                }
                cSelectedProfile aSelection;
                ReadFromFile(aSelection,aSelectedProfile);
                mProfileName = aSelection.mNameProfile;
                CheckProfileName(mProfileName);   // comes from a file, not from the user
            }

            mUserProfile = mDirUserProfile + TheUserProfilePrefix + mProfileName + ".xml";
            const std::string aDefaultProfile = mDirLocalParameters + TheDefaultProfileName;
            mParamProfile = SyncUserProfile(mUserProfile,aDefaultProfile);
        }
        catch (const cProfileError & anError)
        {
            TheProfileFailure = anError.mType + " " + anError.mMes;
            aProfileIsUnusable = ! IsUnreachableProfile(anError.mType);
        }
    }
    if (! TheProfileFailure.empty())
    {
        //  Nothing can be read from or written to a profile, run with the built-in values.
        //  Empty directory marks the profile as unusable, see GlobProfileNames and EditProfile
        //  Warn only when a profile exists but cannot be understood : the user thinks their
        //  settings apply, and they do not. Not being able to reach one stays silent.
        if (aProfileIsUnusable)
           MMVII_USER_WARNING("Invalid user profile, running with the built-in values : " + TheProfileFailure);
        mDirUserProfile = "";
        mUserProfile = "";
        mParamProfile = BuiltInProfile();
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
    //  When there is no usable user configuration directory, no profile file can be
    //  enumerated, only the default one exists, in memory
    if ((! aDirUserProfile.empty()) && IsDirectory(aDirUserProfile))
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

static const tMapProfileValidator & MapProfileValidator()
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


template <typename K, typename V>
std::vector<K> MapKeys(const std::map<K,V>& aMap)
{
    static const std::vector<K> TheKeys = [&aMap]
    {
        std::vector<K> keys;
        for (const auto& [key, _] : aMap)
        {
            keys.push_back(key);
        }
        return keys;
    }();
    return TheKeys;
}

template<typename K, typename V>
std::string MapKeysToS(const std::map<K,V>& aMap)
{
    static const std::string TheKeys = [&aMap]
    {
        std::vector<std::string> keys;
        for (const auto& [key, _] : aMap)
        {
            keys.push_back(ToS(key));
        }
        return ToS(keys);
    }();
    return TheKeys;
}


struct cKeyVal{
    std::string Key;
    std::string Val;
    ARG2007_STRUCT_FIELDS(Key,FieldSem({eTA2007::AllowedValues,MapKeysToS(MapProfileValidator())}),Val);
};


class cAppli_EditProfile : public cMMVII_Appli
{
     public :
        cAppli_EditProfile(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli &);
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;

     private :
         std::string     mCurrent;
         cKeyVal mKeyVal;
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
        << AOpt2007(mCurrent,"SetCurrent","Name of the profile to make current",{{eTA2007::AllowedValues,ToS(GlobProfileNames())}})
        << AOpt2007(mKeyVal,"KeyVal","Set Value to Key: [Key,Value]")
        << AOpt2007(mDelKey,"DelKey","Delete Key",{{eTA2007::AllowedValues,ToS(mParamProfile.Keys())}})
   ;
}



int cAppli_EditProfile::Exe()
{
    //  Startup can go on with default values when the profile is unusable, but this command
    //  was explicitly called to show or modify a profile : here it is an error
    MMVII_INTERNAL_ASSERT_User
    (
        ! mDirUserProfile.empty(),
        eTyUEr::eUnClassedError,
        "Cannot use the user profile : " + TheProfileFailure
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

        for (const auto & [aKey , _] : MapProfileValidator())
        {
            std::string aVal = "";
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
        const auto & aMapValidator = MapProfileValidator();
        const auto aItValidator = aMapValidator.find(mKeyVal.Key);

        if (aItValidator != aMapValidator.end())
        {
            auto& aValidator = aItValidator->second;
            MMVII_INTERNAL_ASSERT_User
            (
                aItValidator->second(mKeyVal.Val),
                eTyUEr::eBadOptParam,
                "Invalid value [" + mKeyVal.Val + "] for profile key [" + mKeyVal.Key
                    + "], expecting " + ( aValidator.Allowed().empty() ? ("type " + aValidator.TypeName()) : aValidator.Allowed() )
            );
            StdOut() << "Setting profile key '" << mKeyVal.Key << "' to value '" << mKeyVal.Val << "'\n";
            aParam.Set(mKeyVal.Key, mKeyVal.Val);
        }
        else
        {
            MMVII_USER_WARNING("Unknown profile key: " + mKeyVal.Key);
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

