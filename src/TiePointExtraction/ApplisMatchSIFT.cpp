#include "MMVII_DeclareCste.h"
#include "MMVII_Tpl_Images.h"

#include <cstdlib>						// C standard library
#include <cstdio>						// C I/O (for sscanf)
#include <cstring>						// string manipulation
#include <fstream>						// file I/O
#include <vector>

#include "ANN/ANN.h"					// ANN declarations
//#include <ANN/ANN.h>					// ANN declarations
#include "KeyPoint.h"
#include "outilsstring/Stringpp.hpp"
#include "AppariementDescripteurs.hpp"
#include "FiltrageVoisinage.h"
#include "FiltrageAppariementUnique.h"
#include "FiltrageSigne.hpp"

#include "ResultMEC.h"
#include "t2d.hpp"
#include "t3d.hpp"

#include "EmprisesClicheCliche.hpp"

void SIFT2RootSIFT(KeyPoint & kp){
  std::vector<double> desc(kp.desc().size());
  //Normalisation L1
  double normL1=0;
  //  std::cout<<"\nnouveau :"<<std::endl;
  //for(KeyPoint::const_iterator i=kp.begin();i!=kp.end();++i){std::cout<<(double)(*i)<<" ";}
  //std::cout<<std::endl;
  /*
  double normL2=0;
  for(KeyPoint::const_iterator i=kp.begin();i!=kp.end();++i){normL2=normL2+std::pow((double)(*i),2.0);}
  normL2=sqrt(normL2);
  //std::cout<<"L2 = "<<normL2<<std::endl;
  double dtmp;
  for(KeyPoint::iterator i=kp.begin();i!=kp.end();++i){
    dtmp=(((double)(*i))*512.0/normL2);
    if(dtmp>255){std::cerr<<(dtmp)<<std::endl;}
    (*i)=(unsigned char)dtmp;
  }
  //for(KeyPoint::const_iterator i=kp.begin();i!=kp.end();++i){std::cout<<(double)(*i)<<" ";}
  */
  for(KeyPoint::const_iterator i=kp.begin();i!=kp.end();++i){normL1=normL1+fabs((double)(*i));}
  //std::cout<<"L1 = "<<normL1<<std::endl;

  //Take the square root of each element in the SIFT vector. Then the vectors are L2-normalized
  double dtmp;
  for(KeyPoint::iterator i=kp.begin();i!=kp.end();++i){
    dtmp=(double)(*i); 
    dtmp=dtmp/normL1; 
    dtmp=sqrt(dtmp);
    dtmp=dtmp*512.0;
    if(dtmp>255){std::cerr<<(dtmp)<<std::endl;}
    (*i)=(unsigned char)dtmp;
  }
}

void SIFT2RootSIFT(std::vector<KeyPoint> & vkp){
  for(std::vector<KeyPoint>::iterator i=vkp.begin();i!=vkp.end();++i){SIFT2RootSIFT((*i));}
}


//#include "getopt.h"
void Aide(std::ostream & fout, const std::string & adresse_exe){
  fout<<"Appariement : "<<std::endl;
  fout<<adresse_exe<<" ";
  fout<<"--kp1 adresse_data ";
  fout<<"--kp2 adresse_query ";
  fout<<"[--result adresse_result ] ";
  fout<<"[--filtrevoisinage (par defaut inactif)] ";
  fout<<"[--filtresigne (par defaut inactif)] ";
  fout<<"[--filtreappariementunique (par defaut inactif)] ";
  fout<<"[--filtredoublons (par defaut inactif)] ";
  fout<<"[--filtreratio (par defaut 0.6) ]";
  fout<<"[--maxptsavisiterkdtree (par defaut 200) ]";
  fout<<"[--pardalle (par defaut 0 = tout) ]";
  fout<<"[--rootsift (par defaut inactif) ]";
  fout<<"[--emprises adresse_emprises] ";
  fout<<"[--verbose] ";
  fout<<"[--sigmamin (pardefaut -1)] ";
  fout<<"[--sigmamax (pardefaut -1)] ";
  fout<<"[--octavemin (pardefaut -2)] ";
  fout<<"[--octavemax (pardefaut -2)] ";
  fout<<"[--nbdecimales (par defaut 4)] ";
  fout<<std::endl;
};


namespace MMVII
{



class cAppli_MatchSIFT : public cMMVII_Appli
{
     public :
        cAppli_MatchSIFT(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec,bool isBasic);
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;
     private :
  std::string mAdresseKpts1, mAdresseKpts2; 
  std::string m_adresse_result;
  bool m_usefiltrevoisinage=false;
  bool m_usefiltresigne=false;
  bool m_usefiltreappariementunique=false;
  bool m_usefiltredoublons=false;
  double m_filtreratio=0.6;
  int m_maxPtsAVisiterDansKDTree=200;
  std::string m_adresse_emprises_predicteur;
  bool m_verbose=false;
  bool m_rootsift=false;
  int m_tailledalle=0;
  double m_sigmamin=-1;
  double m_sigmamax=-1;
  double m_octavemin=-2;
  double m_octavemax=-2;
  int m_nbdecimales=4;
};

cCollecSpecArg2007 & cAppli_MatchSIFT::ArgObl(cCollecSpecArg2007 & anArgObl)
{
   anArgObl
     <<   Arg2007(mAdresseKpts1,"Path to keypoint input file 1",{})
     <<   Arg2007(mAdresseKpts2,"Path to keypoint input file 2",{})
     //<<   Arg2007(mNameIn,"Name of input file",{{eTA2007::MPatFile,"0"},eTA2007::FileImage})
   //          <<   Arg2007(mScale,"Scaling factor",{})
   ;
  return anArgObl;
}

cCollecSpecArg2007 & cAppli_MatchSIFT::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
  
  anArgOpt
    //    << AOpt2007(m_NameOut,"Name of output file","Out",{})

    << AOpt2007(m_adresse_result,"Result","Path to output result file",{})
    << AOpt2007(m_usefiltrevoisinage,"FiltreVoisinage","Use FiltreVoisinage",{})
    << AOpt2007(m_usefiltresigne,"FiltreSigne","Use FiltreSigne",{})
    << AOpt2007(m_usefiltreappariementunique,"FiltreAppariementUnique","Use FiltreAppariementUnique",{})
    << AOpt2007(m_usefiltredoublons,"FilteDoublons","Use Filtredoublons",{})
    << AOpt2007(m_filtreratio,"FiltreRatio","Use Filtreratio",{})
    << AOpt2007(m_maxPtsAVisiterDansKDTree,"maxPtsAVisiterDansKDTree","maxPtsAVisiterDansKDTree",{})
    << AOpt2007(m_tailledalle,"tailledalle","tailledalle",{})
    << AOpt2007(m_rootsift,"RootSift","Use rootsift",{})
    << AOpt2007(m_adresse_emprises_predicteur,"Emprises","Adresse emprise predicteur",{})
    << AOpt2007(m_verbose,"Verbose","Verbose",{})
    << AOpt2007(m_sigmamin,"SigmaMin","Sigma min",{})
    << AOpt2007(m_sigmamax,"SigmaMax","Sigma max",{})
    << AOpt2007(m_octavemin,"OctaveMin","Octave min",{})
    << AOpt2007(m_octavemax,"OctaveMax","Octave max",{})
    << AOpt2007(m_nbdecimales,"NbDecimales","NbDecimales",{})
    ;
  
  

   return anArgOpt;
}



  
int cAppli_MatchSIFT::Exe()
{
  //  std::cout<<"ici0 "<<mNameIn<<std::endl;
  /*
  if (RunMultiSet(0,0))
   {
       return ResultMultiSet();
   }
  */
  // std::cout<<"ici1 "<<mNameIn<<std::endl;
  //std::cout<<"ici1 "<<mNameOut<<std::endl;
  stringpp adresse_result(m_adresse_result);
  
  bool usefiltrevoisinage=m_usefiltrevoisinage;//false;
  bool usefiltresigne=m_usefiltresigne;//false;
  bool usefiltreappariementunique=m_usefiltreappariementunique;//false;
  bool usefiltredoublons=m_usefiltredoublons;//false;
  double filtreratio=m_filtreratio;//0.6;
  int maxPtsAVisiterDansKDTree=m_maxPtsAVisiterDansKDTree;//200;
  stringpp adresse_emprises_predicteur=m_adresse_emprises_predicteur;
  bool verbose=m_verbose;//false;
  bool rootsift=m_rootsift;//false;
  int tailledalle=m_tailledalle;//0;
  double sigmamin=m_sigmamin;//-1;
  double sigmamax=m_sigmamax;//-1;
  double octavemin=m_octavemin;//-2;
  double octavemax=m_octavemax;//-2;
  int nbdecimales=m_nbdecimales;//4;

  bool selectionselonscale=false;
  if(octavemin!=-2){sigmamin=1.6*std::pow(2.0,(double)octavemin);}
  if(octavemax!=-2){sigmamax=1.6*std::pow(2.0,(double)(octavemax+1.0));}
  if(sigmamax>0){selectionselonscale=true;}else{sigmamax=10000000000;}
  if(sigmamin>0){selectionselonscale=true;}
    

  //----------------------------------------
  //Debut des traitements

  if(mAdresseKpts1.empty()){
    if(verbose)std::cout<<"Il faut indiquer l'adresse des points extraits de l'image data..."<<std::endl;
    return !EXIT_SUCCESS;
  }
  if(mAdresseKpts2.empty()){
    if(verbose)std::cout<<"Il faut indiquer l'adresse des points extraits de l'image query..."<<std::endl;
    return !EXIT_SUCCESS;
  }
  
  //--------------------------------------
  //Lecture des points d'interet
  if(verbose)std::cout<<"Lecture des points data..."<<std::flush;
  std::vector<KeyPoint> kpdata,kpquery;
  if(!selectionselonscale){
    if(!OutilsIOKeyPoint::ChargerKeyBin(mAdresseKpts1,kpdata)){return !EXIT_SUCCESS;}
    if(verbose)std::cout<<"OK"<<std::endl<<"Lecture des points query..."<<std::flush;
    if(!OutilsIOKeyPoint::ChargerKeyBin(mAdresseKpts2,kpquery)){return !EXIT_SUCCESS;}
    if(verbose)std::cout<<"OK"<<std::endl<<std::flush;
  }else{
    if(!OutilsIOKeyPoint::ChargerKeyBin(mAdresseKpts1,kpdata,sigmamin,sigmamax)){return !EXIT_SUCCESS;}
    if(verbose)std::cout<<"OK"<<std::endl<<"Lecture des points query..."<<std::flush;
    if(!OutilsIOKeyPoint::ChargerKeyBin(mAdresseKpts2,kpquery,sigmamin,sigmamax)){return !EXIT_SUCCESS;}
    if(verbose)std::cout<<"OK"<<std::endl<<std::flush;
  }

  //--------------------------------------
  //Lecture des emprises (A faire en binaire par la suite ???) et selection des points
  TPolygone2D<double> emprise_query_dans_data,emprise_data_dans_query;
  if(!adresse_emprises_predicteur.empty() && EmprisesClicheCliche::LectureEmprisesTexte(adresse_emprises_predicteur,emprise_query_dans_data,emprise_data_dans_query)){
    //Selection des points data situes dans l'emprise
    if(verbose)std::cout<<"Selection des points data appartenant a l'emprise..."<<std::flush;
    std::vector<KeyPoint> kptmp;
    kptmp.reserve(kpdata.size());
    for(std::vector<KeyPoint>::const_iterator it=kpdata.begin();it!=kpdata.end();++it){
      if(emprise_query_dans_data.contient((*it).x(),(*it).y())) kptmp.push_back(*it);
    }
    if(kptmp.empty()){if(verbose){std::cout<<"Pas assez de points data"<<std::endl;} return !EXIT_SUCCESS;}
    kpdata=kptmp;

    //Selection des points query situes dans l'emprise
    if(verbose)std::cout<<"OK"<<std::endl<<"Selection des points query appartenant a l'emprise..."<<std::flush;
    kptmp.clear();  kptmp.reserve(kpquery.size());
    for(std::vector<KeyPoint>::const_iterator it=kpquery.begin();it!=kpquery.end();++it){
      if(emprise_data_dans_query.contient((*it).x(),(*it).y())) kptmp.push_back(*it);
    }
    if(kptmp.empty()){if(verbose){std::cout<<"Pas assez de points query"<<std::endl;} return !EXIT_SUCCESS;}
    kpquery=kptmp;

    if(verbose)std::cout<<"OK"<<std::endl;
  }else{if(verbose)std::cout<<"Pas d'emprise..."<<std::endl;}

  //rootsift
  if(rootsift){
    if(verbose){std::cout<<"Preparation Rootsift."<<std::endl;}
    SIFT2RootSIFT(kpdata);
    SIFT2RootSIFT(kpquery);
  }

  ResultMEC<KeyPoint> mec;

  if(tailledalle<1){
    //------------------------------------------------
    //Appariement
    AppariementDescripteurs appdesc;
    appdesc.verbose()=verbose;
    appdesc.ratiodistances()=filtreratio;//ajout 2017
    appdesc.maxPtsAvisiter()=maxPtsAVisiterDansKDTree;//Ajout 2020
    std::vector<unsigned int> iddata_app, idquery_app;
    appdesc.apparie(kpdata,kpquery,iddata_app,idquery_app);
  
    //------------------------------------------------
    //On recupere uniquement les points 2D
    mec.pts1().resize(iddata_app.size());
    mec.pts2().resize(idquery_app.size());
    for(unsigned int n=0;n<mec.size();n++){
      mec.pt1(n)=kpdata[iddata_app[n]].sansdesc();//.xy();
      mec.pt2(n)=kpquery[idquery_app[n]].sansdesc();//.xy();
    }
    if(verbose)std::cout<<mec.card()<<" points apparies."<<std::endl;
  }else{
    //Il faut commencer par calculer l'emprise des points data...
    double Xmin=kpdata.front().x();
    double Ymin=kpdata.front().y();
    double Xmax=Xmin;
    double Ymax=Ymin;
    for(std::vector<KeyPoint>::const_iterator it=kpdata.begin();it!=kpdata.end();++it){
      if((*it).X()<Xmin){Xmin=(*it).X();}else if((*it).X()>Xmax){Xmax=(*it).X();}
      if((*it).Y()<Ymin){Ymin=(*it).Y();}else if((*it).Y()>Ymax){Ymax=(*it).Y();}
    }
    if(verbose){std::cout<<"("<<Xmin<<";"<<Ymin<<") -> ("<<Xmax<<";"<<Ymax<<")"<<std::endl;}
    //if(emprise_query_dans_data.contient((*it).x(),(*it).y())) kptmp.push_back(*it);
    double dX=Xmax-Xmin;
    double dY=Ymax-Ymin;
    if(verbose){std::cout<<dX/tailledalle<<" x "<<dY/tailledalle<<" dalles"<<std::endl;}
    for(double X0=Xmin;X0<=Xmax;X0=X0+tailledalle){
      double X1=X0+tailledalle;
      for(double Y0=Ymin;Y0<=Ymax;Y0=Y0+tailledalle){
	double Y1=Y0+tailledalle;
	if(verbose){std::cout<<" Dalle : ("<<X0<<";"<<Y0<<") -> ("<<X1<<";"<<Y1<<")"<<std::endl;}
	std::vector<KeyPoint> kpdatatmp;
	for(std::vector<KeyPoint>::const_iterator it=kpdata.begin();it!=kpdata.end();++it){
	  if((*it).X()<X0){continue;}
	  if((*it).X()>X1){continue;}
	  if((*it).Y()<Y0){continue;}
	  if((*it).Y()>Y1){continue;}
	  kpdatatmp.push_back((*it));
	}

	if(verbose){std::cout<<" "<<kpdatatmp.size()<<" points data"<<std::endl;}
	AppariementDescripteurs appdesc;
	appdesc.verbose()=verbose;
	appdesc.ratiodistances()=filtreratio;//ajout 2017
	appdesc.maxPtsAvisiter()=maxPtsAVisiterDansKDTree;//Ajout 2020
	std::vector<unsigned int> iddata_app, idquery_app;
	appdesc.apparie(kpdatatmp,kpquery,iddata_app,idquery_app);
	
	//On recupere uniquement les points 2D
	for(unsigned int n=0;n<idquery_app.size();n++){
	  mec.pts1().push_back(kpdatatmp[iddata_app[n]].sansdesc());//.xy();
	  mec.pts2().push_back(kpquery[idquery_app[n]].sansdesc());//.xy();
	}
	if(verbose)std::cout<<" "<<idquery_app.size()<<" points apparies."<<std::endl;

      }
    }
    if(verbose)std::cout<<mec.card()<<" points apparies."<<std::endl;
  }
  

  //------------------------------------------------
  //FiltrageSigne (au cas ou)
  if(usefiltresigne){
    if(verbose)std::cout<<"Filtre signe :"<<std::endl;
    FiltreSigne filtresigne;
    filtresigne.filtre(mec.pts1(),mec.pts2());
    if(verbose)std::cout<<mec.card()<<" points apparies."<<std::endl;
  }
  
  //------------------------------------------------
  //FiltrageVoisinage (au cas ou)
  if(usefiltrevoisinage){
    if(verbose)std::cout<<"Filtre voisinage :"<<std::endl;
    FiltreVoisinage filtrevoisinage;
    filtrevoisinage.filtre(mec.pts1(),mec.pts2()); 
    filtrevoisinage.filtre(mec.pts1(),mec.pts2());
    if(verbose)std::cout<<mec.card()<<" points apparies."<<std::endl;
  }

  //------------------------------------------------
  //FiltrageAppariementUnique (au cas ou)
  if(usefiltreappariementunique){
    if(verbose)std::cout<<"Filtre appariement unique :"<<std::endl;
    FiltreAppariementUnique filtreappariementunique(true);
    filtreappariementunique.filtre(mec.pts1(),mec.pts2()); 
    if(verbose)std::cout<<mec.card()<<" points apparies."<<std::endl;
  }else if(usefiltredoublons){//simple filtrage des doublons dans ce cas...
    FiltreAppariementUnique::FiltrageDoublons(mec.pts1(),mec.pts2());
    if(verbose)std::cout<<mec.card()<<" points apparies."<<std::endl;
  }

  //------------------------------------------------
  //Sauvegarde du resultat
  if(adresse_result.empty()){adresse_result.CombineNoms(mAdresseKpts1,mAdresseKpts2,".result");}
  mec.EcrireTexte(adresse_result,nbdecimales);

  
   return EXIT_SUCCESS;
}

  
  cAppli_MatchSIFT:: cAppli_MatchSIFT(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec,bool isBasic):
    cMMVII_Appli (aVArgs,aSpec),
    m_adresse_result(""),
    m_usefiltrevoisinage(false),
    m_usefiltresigne(false),
    m_usefiltreappariementunique(false),
    m_usefiltredoublons(false),
    m_filtreratio(0.6),
    m_maxPtsAVisiterDansKDTree(200),
    m_adresse_emprises_predicteur(""),
    m_verbose(false),
    m_rootsift(false),
    m_tailledalle(0),
    m_sigmamin(-1),
    m_sigmamax(-1),
    m_octavemin(-2),
    m_octavemax(-2),
    m_nbdecimales(4)

{
}

  //====================  Appli MatchSIFT SPEC =======================================================
  
tMMVII_UnikPApli Alloc_MatchSIFT(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppli_MatchSIFT(aVArgs,aSpec,true));
}

  //Descripteur pour faire connaitre la classe au reste de 2007
cSpecMMVII_Appli  TheSpecMatchSIFT
(
     "MatchSIFT",
     Alloc_MatchSIFT,
     "Match SIFT keypoints",
     {eApF::TieP},
     {eApDT::Image},
     {eApDT::TieP},
      __FILE__
);



};
