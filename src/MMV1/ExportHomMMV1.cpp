#include "V1VII.h"

#include "MMVII_util.h"
#include "MMVII_MeasuresIm.h"
#include "MMVII_Stringifier.h"

namespace MMVII
{

#if (MMVII_KEEP_LIBRARY_MMV1)
cHomogCpleIm  ToMMVII(const cNupletPtsHomologues &  aNUp)
{

        return  cHomogCpleIm(ToMMVII(aNUp.P1()),ToMMVII(aNUp.P2()));
}

/*   ************************************************* */
/*                                                     */
/*         cConvertHomV1                               */
/*                                                     */
/*   ************************************************* */

class cImportHomV1 : public cInterfImportHom
{
      public :
          cImportHomV1(const std::string & aDir,const std::string & aSubDir,const std::string & aExt="dat");

          void GetHom(cSetHomogCpleIm &,const std::string & aNameIm1,const std::string & aNameIm2) const override;
          bool   HasHom(const std::string & aNameIm1,const std::string & aNameIm2) const override;


      private :

          std::string NameHom(const std::string & aNameIm1,const std::string & aNameIm2) const;
          std::string  mDir;
          std::string  mSubDir;
          std::string  mExt;
          std::string  mKHIn;
          cElemAppliSetFile                 mEASF;
          cInterfChantierNameManipulateur * mICNM ;
};


cImportHomV1::cImportHomV1(const std::string & aDir,const std::string & aSubDir,const std::string & aExt) :
   mDir     (aDir) ,
   mSubDir  (aSubDir) ,
   mExt     (aExt),
   mKHIn    (std::string("NKS-Assoc-CplIm2Hom@") + mSubDir  +  std::string("@") +  mExt),
   mEASF    (mDir),
   mICNM    (mEASF.mICNM)
{
}


std::string cImportHomV1::NameHom(const std::string & aNameIm1,const std::string & aNameIm2) const
{
    return mDir +  mICNM->Assoc1To2(mKHIn,aNameIm1,aNameIm2,true);
}

bool   cImportHomV1::HasHom(const std::string & aNameIm1,const std::string & aNameIm2) const
{
    return ExistFile(NameHom(aNameIm1,aNameIm2));
}

void  cImportHomV1::GetHom(cSetHomogCpleIm & aPackV2,const std::string & aNameIm1,const std::string & aNameIm2) const
{
     aPackV2.Clear();
     ElPackHomologue aPackV1 = ElPackHomologue::FromFile(NameHom(aNameIm1,aNameIm2));

     for (ElPackHomologue::tIter aItV1 = aPackV1.begin() ; aItV1!=aPackV1.end() ; aItV1++)
     {
         aPackV2.Add(ToMMVII(*aItV1));
     }
}

/*   ************************************************* */
/*                                                     */
/*                 cInterfImportHom                    */
/*                                                     */
/*   ************************************************* */

cInterfImportHom * cInterfImportHom::CreateImportV1(const std::string&aDir,const std::string&aSubDir,const std::string&aExt)
{
        return new cImportHomV1(aDir,aSubDir,aExt);
}
#else // MMVII_KEEP_LIBRARY_MMV1

/**
 * @brief The cMMVII_Implem_V1ImportHom class
 *
 * A class for importing V1-Homologous in MMVII w/o using the V1 library
 */
class cMMVII_Implem_V1ImportHom : public cInterfImportHom
{
      public :
          cMMVII_Implem_V1ImportHom(const std::string & aDir,const std::string & aSubDir,const std::string & aExt="txt");

          void GetHom(cSetHomogCpleIm &,const std::string & aNameIm1,const std::string & aNameIm2) const override;
          bool   HasHom(const std::string & aNameIm1,const std::string & aNameIm2) const override;


      private :

          std::string NameHom(const std::string & aNameIm1,const std::string & aNameIm2) const;
          std::string  mDir;
          std::string  mSubDir;
          std::string  mExt;
         // std::string  mKHIn;
};

cMMVII_Implem_V1ImportHom::cMMVII_Implem_V1ImportHom(const std::string & aDir,const std::string & aSubDir,const std::string & aExt) :
    mDir    (aDir),
    mSubDir (aSubDir),
    mExt    (aExt)
{
}


std::string cMMVII_Implem_V1ImportHom::NameHom(const std::string & aNameIm1,const std::string & aNameIm2) const
{
    return mDir
            + std::string("Homol") + mSubDir + StringDirSeparator()
            + "Pastis" + aNameIm1 +  StringDirSeparator()
            + aNameIm2 + ".txt"
    ;
}

bool   cMMVII_Implem_V1ImportHom::HasHom(const std::string & aNameIm1,const std::string & aNameIm2) const
{
    bool hasHom =  ExistFile(NameHom(aNameIm1,aNameIm2));

    return hasHom;
}

void cMMVII_Implem_V1ImportHom::GetHom(cSetHomogCpleIm & aSetH,const std::string & aNameIm1,const std::string & aNameIm2) const
{
   std::string aNameFile = NameHom(aNameIm1,aNameIm2);
   std::ifstream aInputFil(aNameFile);

   std::string aLine;
   int aCptLine=0;
   while (std::getline(aInputFil, aLine))
   {
       std::istringstream iss(aLine);
       std::vector<double> aVCoord;

       for (int aK=0 ; aK<4 ; aK++)
       {
           tREAL8 aVal = GetV<tREAL8>(iss,aNameFile,aCptLine);
           aVCoord.push_back(aVal);

       }
       cPt2dr aPIm1(aVCoord[0],aVCoord[1]);
       cPt2dr aPIm2(aVCoord[2],aVCoord[3]);

       aSetH.Add(cHomogCpleIm(aPIm1,aPIm2));
       aCptLine++;
   }

   // aSetH.InitFromFile(NameHom(aNameIm1,aNameIm2));
   //
   // StdOut() << "GET HOM " << aNameIm1 << "/" << aNameIm2 << " NBH=" << aCptLine << "\n";
}


cInterfImportHom * cInterfImportHom::CreateImportV1(const std::string&aDir,const std::string&aSubDir,const std::string&aExt)
{
        MMVII_INTERNAL_ASSERT_User_UndefE(aExt=="txt","CreateImportV1 w/o V1-Lib requires V1-format");
        return new cMMVII_Implem_V1ImportHom(aDir,aSubDir,aExt);
}

#endif // MMVII_KEEP_LIBRARY_MMV1


cInterfImportHom::~cInterfImportHom()
{
}



};
