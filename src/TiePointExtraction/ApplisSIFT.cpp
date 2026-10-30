#include "MMVII_DeclareCste.h"
#include "MMVII_Tpl_Images.h"

#include<string>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<algorithm>

#include "sift-driver.hpp"

namespace MMVII
{



class cAppli_SIFT : public cMMVII_Appli
{
     public :
        cAppli_SIFT(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec,bool isBasic);
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;
     private :
        ParametrageSiftDriver mParametres;
        std::string mNameIn;  ///< Input image name
        std::string mNameOut; ///< Output file name
        std::string mOutputFilenamePrefix, mInputKeyPointsToBeDescribed;

};

cCollecSpecArg2007 & cAppli_SIFT::ArgObl(cCollecSpecArg2007 & anArgObl)
{
   anArgObl
     <<   Arg2007(mNameIn,"Name of input file",{})
     //<<   Arg2007(mNameIn,"Name of input file",{{eTA2007::MPatFile,"0"},eTA2007::FileImage})
   //          <<   Arg2007(mScale,"Scaling factor",{})
   ;
  return anArgObl;
}

cCollecSpecArg2007 & cAppli_SIFT::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
  anArgOpt
    << AOpt2007(mNameOut,"Out","Name of output file",{eTA2007::HDV})
    << AOpt2007(mParametres.verbose(),"Verbose","Verbosity",{eTA2007::HDV})
    << AOpt2007(mParametres.sift().firstoctave(),"FirstOctave","First octave",{eTA2007::HDV})
    << AOpt2007(mParametres.sift().nboctaves(),"Octaves","Nb octaves",{eTA2007::HDV})
    << AOpt2007(mParametres.sift().nblevels(),"Levels","Nb levels",{eTA2007::HDV})
    << AOpt2007(mParametres.sift().threshold(),"Threshold","Threshold",{eTA2007::HDV})
    << AOpt2007(mParametres.sift().edgeThreshold(),"EdgeThreshold","Edge-threshold",{eTA2007::HDV})
    //    << AOpt2007(mParametres.pretraitement().canal(),"Channel","Channel",{})
    << AOpt2007(mParametres.pretraitement().maxim(),"Max","Valeur maximale 255 ou autre (par defaut -1 = vraie valeur maximale)",{eTA2007::HDV})
    << AOpt2007(mParametres.pretraitement().minim(),"Min","Valeur minimale 255 ou autre (par defaut -1 = vraie valeur maximale)",{eTA2007::HDV})
    << AOpt2007(mParametres.sift().signe(),"Sign","Sign  : 0 > et < ; -1 < ; 1 >",{eTA2007::HDV})
    << AOpt2007(mOutputFilenamePrefix,"OutPrefix","OutputFileNamePrefix",{eTA2007::HDV})
    << AOpt2007(mInputKeyPointsToBeDescribed,"InPts","Input keypoints to be described",{eTA2007::HDV})
    << AOpt2007(mParametres.sift().siftdense(),"DenseStep","Dense SIFT spatial step",{eTA2007::HDV})
    << AOpt2007(mParametres.sift().siftdensescalefixe(),"DenseScale","Dense SIFT spatial step",{eTA2007::HDV})
    << AOpt2007(mParametres.sift().binary(),"Binary","Binary format",{eTA2007::HDV})
    << AOpt2007(mParametres.pretraitement().usegradient(),"UseGradient","Use gradient image",{eTA2007::HDV})
    << AOpt2007(mParametres.sift().demidescripteur(),"HalfDesc","Half descriptors",{eTA2007::HDV})
    << AOpt2007(mParametres.dallage().taille_souhaitee_dalle(),"Tile","Tile size",{eTA2007::HDV})
    << AOpt2007(mParametres.sift().nodescr(),"NoDescr","No descriptor",{eTA2007::HDV})
    << AOpt2007(mParametres.sift().noorient(),"NoOrient","No orientation",{eTA2007::HDV})
  ;

    /*     anArgOpt
             << AOpt2007(mDilate,"Dilate","Dilatation od gaussian",{eTA2007::HDV})
             << AOpt2007(mForceGray,"FG","Force gray image",{eTA2007::HDV})
             << AOpt2007(mPrefOut,"PrefOut","Prefix for output image",{eTA2007::HDV})
             << AOpt2007(mNameOut,CurOP_Out,"Name of out file, Def=\"$PrefOut\"+Input",{})
   ;
  */
  

   return anArgOpt;
}


  void DecouverteMMVII(const std::string NameIm){
      cDataFileIm2D aFileIn= cDataFileIm2D::Create(NameIm,eForceGray::No);
  std::cout<<aFileIn.Sz()<<std::endl;
  std::cout<<aFileIn.Sz().x()<<std::endl;
  std::cout<<aFileIn.Sz().y()<<std::endl;

  //Apparemment il faut definir la taille de l'image puis decider a partir d'ou on lit...
  cIm2D<tREAL4> aImIn(aFileIn.Sz());
  aImIn.Read(aFileIn,cPt2di(0,0));
  //  aImIn.Read(aFileIn,cBox2di(cPt2di(0,0),cPt2di(10,5)));

  // aImIn.DIm().ExtractRawData2D();
  //(y*longueurligne+x)*taillemem

  /*
  
   for(int i=0;i<aImIn.DIm().Sz().x();i++){
     for(int j=0;j<aImIn.DIm().Sz().y();j++){
       std::cout<<aImIn.DIm().GetV(cPt2di(i,j))<<" ";
     }
     std::cout<<std::endl;
   }
   std::cout<<"AAA"<<std::endl;
   tREAL4** pt=aImIn.DIm().ExtractRawData2D();
   for(int i=0;i<aImIn.DIm().Sz().x();i++){
     for(int j=0;j<aImIn.DIm().Sz().y();j++){
       std::cout<<pt[j][i]<<" ";
     }
     std::cout<<std::endl;
   }


   aImIn.DIm().ToFile("aaaa.tif");
   cIm2D<tREAL4> aImIn1(cPt2di(100,100));
   aImIn1.Read(aFileIn,cPt2di(10,10));
   aImIn1.DIm().ToFile("aaaa.tif");
   pt=aImIn1.DIm().ExtractRawData2D();
   for(int c=0;c<aImIn1.DIm().Sz().x();c++){
     for(int l=0;l<aImIn.DIm().Sz().y();l++){
       std::cout<<pt[l][c]<<" ";
       }
     std::cout<<std::endl;
   }
   

   for(int c=0;c<aImIn1.DIm().Sz().x();c++){
     for(int l=0;l<aImIn1.DIm().Sz().y();l++){
       std::cout<<aImIn1.DIm().GetV(cPt2di(c,l))<<"="<<pt[l][c]<<" ";
     }
     std::cout<<std::endl;
   }
*/

  }

  
int cAppli_SIFT::Exe()
{
  //  std::cout<<"ici0 "<<mNameIn<<std::endl;
  /*
  if (RunMultiSet(0,0))
   {
       return ResultMultiSet();
   }
  */
  //  std::cout<<"ici1 "<<mNameIn<<std::endl;
  //  std::cout<<"ici1 "<<mNameOut<<std::endl;
  


  stringpp outputFilename(mNameOut);

  
  //Gestion des options / arguments : verifications que 
  if(mParametres.sift().nboctaves()<1 && mParametres.sift().nboctaves()!=-1){std::cerr<<"Number of octaves cannot be smaller than one."<<std::endl; return !EXIT_SUCCESS;}
  if(mParametres.sift().nblevels()< 1) {std::cerr<<"Number of levels cannot be smaller than one."<<std::endl; return !EXIT_SUCCESS;}
  //  if(abs(mParametres.sift().signe())>2){mParametres.sift().signe()=0;}
  if(!mInputKeyPointsToBeDescribed.empty()){mParametres.sift().lireKeyPointsADecrire(mInputKeyPointsToBeDescribed);}
  if(outputFilename.size() !=0 && mOutputFilenamePrefix.size() !=0) {
    std::cerr<<"OutputFileName OR OutputFilenamePrefix, but not both"<<std::endl; return !EXIT_SUCCESS;
  }

  //Au cas ou l'on souhaite faire du sift dense, on peut revenir sur certains parametres...
  if(mParametres.sift().issiftdense()){
    mParametres.sift().update_cas_siftdense();
    if(mParametres.verbose()){std::cout<<"siftdense pour un pas de "<<mParametres.sift().siftdense()<<std::endl;}
    if(mParametres.verbose()){std::cout<<"calcul ne concernant que l'octave "<<mParametres.sift().firstoctave()<<std::endl;}
  }

  
  stringpp adresse_image(mNameIn);
  stringpp outputFilenamePrefix;
  if(outputFilename.empty()){
    outputFilename=adresse_image ;
    if(!outputFilenamePrefix.empty()){
      // if we specify an output directory, then extract the basename
      outputFilename.adressefromDirFile(outputFilenamePrefix.Dir(),adresse_image.File());
    }
    // remove image extension, add .key
    outputFilename.SetExtension(".key");
  }
  
  //------------Debut du traitement d'une image----------------
  //Creation du flux de sortie
  std::ofstream fkey(outputFilename.c_str(),std::ios::out | std::ios::binary);
  if(!fkey.good()){std::cout<<"Probleme en ecriture du fichier "<<outputFilename<<std::endl; return 1;}
  OutilsIOKeyPoint::FluxSortieKeypoint_fichierbinaire ffkey(fkey);

   
  //Traitement sift de l'image
  SiftDriver siftdriver(mParametres);
  siftdriver.traiter_image(adresse_image,ffkey);
  
  //On peut maintenant ecrire le nombre de points dans le fichier
  ffkey.finir();
  std::cout<<ffkey.size()<<" points"<<std::endl;
  fkey.close();
  //------------Fin du traitement d'une image------------------

    /*
    { "verbose",            no_argument,            NULL,              'v' },
    { "help",               no_argument,            NULL,              'h' },
    { "output",             required_argument,      NULL,              'o' },
    { "prefix",             required_argument,      NULL,              'p' },
    { "first-octave",       required_argument,      NULL,              'f' },
    { "keypointsadecrire",  required_argument,      NULL,              'k' },
    { "siftdense",          required_argument,      NULL,              'D' },
    { "siftdensescalefixe", required_argument,      NULL,              'd' },
    { "octaves",            required_argument,      NULL,              'O' },
    { "levels",             required_argument,      NULL,              'S' },
    { "threshold",          required_argument,      NULL,              't' },
    { "edge-threshold",     required_argument,      NULL,              'e' },
    { "signe",              required_argument,      NULL,              's' },//Ajout signe < >
    { "canal",              required_argument,      NULL,              'c' },//Ajout canal
    { "gradient",           no_argument,            NULL,              'g' },//utilisation du gradient
    { "demidescripteur",    no_argument,            NULL,              'r' },//demi-descripteur
    { "maxim",              required_argument,      NULL,              'M' },//Valeur maximale 255 ou autre (par defaut -1 = vraie valeur maximale) ?
    { "minim",              required_argument,      NULL,              'm' },//Valeur maximale 255 ou autre (par defaut -1 = vraie valeur maximale) ?
    { "dallage",            required_argument,      NULL,              'T' },
    //    { "binary",          no_argument,            NULL,              'b' }, 
    { "no-descriptors",     no_argument,            &(parametres.sift().nodescr()),          1   },
    { "no-orientations",    no_argument,            &(parametres.sift().noorient()),         1   },
    //    { "save-gss",        no_argument,            &parametres.savegss(),          1   },
    { NULL,              0,                      NULL,              0   }
    */
  
   return EXIT_SUCCESS;
}

  
  cAppli_SIFT:: cAppli_SIFT(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec,bool isBasic):
    cMMVII_Appli (aVArgs,aSpec),mParametres(),mNameOut(""),mOutputFilenamePrefix(""),mInputKeyPointsToBeDescribed("")
{
}

  //====================  Appli SIFT SPEC =======================================================
  
tMMVII_UnikPApli Alloc_SIFT(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppli_SIFT(aVArgs,aSpec,true));
}

  //Descripteur pour faire connaitre la classe au reste de 2007
cSpecMMVII_Appli  TheSpecSIFT
(
     "SIFT",
     Alloc_SIFT,
     "Detect SIFT keypoints",
     {eApF::TieP},
     {eApDT::Image},
     {eApDT::TieP},
      __FILE__
);



};

