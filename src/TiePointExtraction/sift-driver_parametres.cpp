#include "sift-driver_parametres.hpp"
#include "outilsstring/Stringpp.hpp"
#include "KeyPoint.h"

const bool ParametrageSiftDriver::LectureFichierParametres(std::string & adresse_fichier_parametres){
  std::ifstream ficparam(adresse_fichier_parametres.c_str(),std::ios::in);
  if(!LectureFichierParametres(ficparam)){
    std::cout<<"Probleme en lecture du fichier "<< adresse_fichier_parametres<<std::endl; 
    return false;
  }
  return true;
}

const bool ParametrageSiftDriver::LectureFichierParametres(std::ifstream & ficparam){
  if(!ficparam.good()){return false;}	
  while(true){
    stringpp stmp;
    ficparam>>stmp;
    if(ficparam.eof()){break;}
    if(stmp.empty()){continue;}
    stmp.ToUpCase();
    stmp.supprimeavantapres_ep();

    //parametrage pretraitement
    if(stmp=="CANAL"){ficparam>>stmp; _param_pretraitement.canal()=stmp.atoi();}
    else if(stmp=="USEGRADIENT"){_param_pretraitement.usegradient()=true;}
    else if(stmp=="MAXIM"){ficparam>>stmp; _param_pretraitement.maxim()=stmp.atof();}
    else if(stmp=="MINIM"){ficparam>>stmp; _param_pretraitement.minim()=stmp.atof();}

    //parametrage dallage
    else if(stmp=="TAILLESOUHAITEEDALLE"){ficparam>>stmp; _param_dallage.taille_souhaitee_dalle()=stmp.atoi();}
    else if(stmp=="TAILLEMINIMALEDALLE"){ficparam>>stmp; _param_dallage.taille_minimale_dalle()=stmp.atoi();}
    else if(stmp=="MARGE"){ficparam>>stmp; _param_dallage.marge()=stmp.atoi();}

    //parametrage sift
    else if(stmp=="FIRSTOCTAVE"){ficparam>>stmp; _param_sift.firstoctave()=stmp.atoi();}
    else if(stmp=="NBOCTAVES"){ficparam>>stmp; _param_sift.nboctaves()=stmp.atoi();}
    else if(stmp=="NBLEVELS"){ficparam>>stmp; _param_sift.nblevels()=stmp.atoi();}
    else if(stmp=="THRESHOLD"){ficparam>>stmp; _param_sift.threshold()=stmp.atof();}
    else if(stmp=="EDGETHRESHOLD"){ficparam>>stmp; _param_sift.edgeThreshold()=stmp.atof();}
    else if(stmp=="SIGNE"){ficparam>>stmp; _param_sift.signe()=stmp.atoi();}   
    else if(stmp=="NODESCR"){_param_sift.nodescr()=true;}   
    else if(stmp=="NOORIENT"){_param_sift.noorient()=true;}   
    else if(stmp=="TAILLEGRILLE_LIGNES"){ficparam>>stmp; }   
 
   else if(stmp=="VERBOSE"){verbose()=true;}
  }

  ficparam.close();
  return true;

}


const bool ParametrageSift::lireKeyPointsADecrire(const std::string & keypointsFilename){
  stringpp ext=stringpp(keypointsFilename).GetExtension().UpCase();
  if(ext==stringpp("TXT")){
    return ParametrageSift::lireKeyPointsADecrire_formattxt(keypointsFilename);
  }else if(ext==stringpp("key")){
    return ParametrageSift::lireKeyPointsADecrire_formatbin(keypointsFilename);
  }
  return ParametrageSift::lireKeyPointsADecrire_formattxt(keypointsFilename);
}

const bool ParametrageSift::lireKeyPointsADecrire_formatbin(const std::string & keypointsFilename){
  std::ifstream fkp(keypointsFilename.c_str(),std::ios::in | std::ios::binary);
  if(!fkp.good()){
    std::cout<<"Impossible d'acceder au fichier '"<<keypointsFilename<<"'"<<std::endl;
    return false;
  }
  return lireKeyPointsADecrire_formatbin(fkp);
}


const bool ParametrageSift::lireKeyPointsADecrire_formattxt(const std::string & keypointsFilename){
  std::ifstream fkp(keypointsFilename.c_str(),std::ios::in);
  if(!fkp.good()){
    std::cout<<"Impossible d'acceder au fichier '"<<keypointsFilename<<"'"<<std::endl;
    return false;
  }
  return lireKeyPointsADecrire_formattxt(fkp);
}

const bool ParametrageSift::lireKeyPointsADecrire_formatbin(std::ifstream & fkp){
  std::vector<KeyPoint> vkp;
  bool ok=OutilsIOKeyPoint::ChargerKeyBin(fkp,vkp);
  _KeyPointsADecrire.resize(vkp.size());
  std::vector<T3D<float> >::iterator i=_KeyPointsADecrire.begin();
  for(std::vector<KeyPoint>::const_iterator j=vkp.begin();j!=vkp.end();++j,++i){
    (*i).x()=(*j).x();
    (*i).y()=(*j).y();
    (*i).z()=(*j).sigma();
  }
  return ok;
}

const bool ParametrageSift::lireKeyPointsADecrire_formattxt(std::ifstream & fkp){
  unsigned int nbpts, dim;//dim=2 (on ne fournit que x et y) ou 3 (on fournit (x,y,sigma))
  fkp>>nbpts>>dim;
  _KeyPointsADecrire.resize(nbpts);
  if(dim==3){
    for(std::vector<T3D<float> >::iterator i=_KeyPointsADecrire.begin();i!=_KeyPointsADecrire.end();++i){
      fkp>>(*i).x()>>(*i).y()>>(*i).z();
    }
  }else if(dim==2){
    for(std::vector<T3D<float> >::iterator i=_KeyPointsADecrire.begin();i!=_KeyPointsADecrire.end();++i){
      fkp>>(*i).x()>>(*i).y();
      (*i).z()=2.6;//correspond a l'octave 0...
    }
  }else{return false;}
  return true;
}

