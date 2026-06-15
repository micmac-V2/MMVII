#include "cMMVII_Appli.h"
#include <clocale>
#include "MMVII_Sys.h"

using namespace MMVII;


#ifdef MMVII_KEEP_MMV1_IMAGE
static constexpr const char* ENV_MMVII_USE_MMV1_IMAGE = "MMVII_USE_MMV1_IMAGE";
extern bool mmvii_use_mmv1_image;
#endif


static bool containsCaseInsensitive(const std::string& str,
                             const std::string& search)
{
    return std::search(
               str.begin(), str.end(),
               search.begin(), search.end(),
               [](char c1, char c2)
               {
                   return std::tolower(static_cast<unsigned char>(c1))
                   == std::tolower(static_cast<unsigned char>(c2));
               }
               ) != str.end();
}

static bool cmd_match(const std::string& aCmd, const std::string& (cSpecMMVII_Appli::*method)() const, const std::string &descr)
{
    const auto &allSpecs = cSpecMMVII_Appli::VecAll();
    std::vector<const cSpecMMVII_Appli*> matchName;
    std::copy_if(
        allSpecs.begin(), allSpecs.end(),
        std::back_inserter(matchName),
        [&](const cSpecMMVII_Appli* aSpec)
        {
            return containsCaseInsensitive((aSpec->*method)(), aCmd);
        }
    );
    if (matchName.size()) {
        StdOut() << "** Command '" << aCmd << "' not found, but the following commands contain '" << aCmd << "' in their " << descr << ":\n";
        for (const auto& aSpec : matchName) {
            StdOut()  << aSpec->Name() << " => " << aSpec->Comment() << std::endl;
        }
        return true;
    }
    return false;
}



int main(int argc, char ** argv)
{
#ifdef MMVII_KEEP_MMV1_IMAGE
    char *env_mmv1_image = getenv(ENV_MMVII_USE_MMV1_IMAGE);
    mmvii_use_mmv1_image = false;
    if (env_mmv1_image != nullptr) {
        if (UCaseEqual(env_mmv1_image,"on") || UCaseEqual(env_mmv1_image,"true") || UCaseEqual(env_mmv1_image,"1"))
            mmvii_use_mmv1_image = true;
    }
#endif

   std::setlocale(LC_ALL, "C");
   // std::setlocale(LC_ALL, "en_US.UTF-8");

   cMMVII_Appli::InitMMVIIDirs(MMVII_CanonicalRootDirFromExec());
    
   if (argc>1)
   {
      std::string aNameCom = argv[1];

      // Recherche la specif correspondant au nom de commande
      const cSpecMMVII_Appli*  aSpec = cSpecMMVII_Appli::SpecOfName(aNameCom,true);

      // Execute si match
      if (aSpec)
      {
         std::vector<std::string> aVArgs;
         for (int aK=0 ; aK<argc; aK++)
             aVArgs.push_back(argv[aK]);
         int aRes =  aSpec->AllocExecuteDestruct(aVArgs);
         return aRes;
      }

      // Command not found
      if (! UCaseEqual(aNameCom,"help") && aNameCom !=  "?")
      {
          if (aNameCom.size() && aNameCom.back() == '?')
          {
              aNameCom.pop_back();
          } else {
              if (cmd_match(aNameCom, &cSpecMMVII_Appli::Name, "name"))
                  return EXIT_SUCCESS;
          }
          if (cmd_match(aNameCom, &cSpecMMVII_Appli::Comment, "comment"))
              return EXIT_SUCCESS;
      }
       StdOut() << "Command '" << aNameCom << "' not found.\n";
       StdOut() << "Use 'MMVII' without any argument to get the list of available commands.\n";
       return EXIT_FAILURE;
   }

   // Affiche toutes les commandes
   for (const auto & aSpec : cSpecMMVII_Appli::VecAll())
   {
      StdOut()  << aSpec->Name() << " => " << aSpec->Comment() << std::endl;
   }

#ifdef MMVII_KEEP_MMV1_IMAGE
   StdOut()
       << "\n"
       << " >>> MMVII is using " << (mmvii_use_mmv1_image ? "<< MicMac v1 >>" : "<< GDal >>") << " to read/write image file.\n"
       << " >>> (Env var '" << ENV_MMVII_USE_MMV1_IMAGE << "' is " << (mmvii_use_mmv1_image ? "" : "NOT ") << "set to 'on', 'true' or '1')" << std::endl;
#endif

   return EXIT_SUCCESS;
}



