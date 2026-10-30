#ifndef _SIFTDRIVER_PARAMETRES_HPP_
#define _SIFTDRIVER_PARAMETRES_HPP_

#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>

//=========================================================
//Classes de gestion des differents parametres utilises lors de la mise en oeuvre de sift
//  - ParametrageDallage : parametres lies au dallage des images
//  - ParametragePretraitement : parametres lies a divers pretraitements (choix du canal, gradient ? , reetalement)
//  - ParametrageSift : parametres lies directement a sift
//  - ParametrageSiftDriver : regroupement des parametres precedents
//=========================================================

#include "Parametres_Pretraitement_Dallage.hpp"
#include "t3d.hpp"



//=========================================================
//ParametrageSift
//
//Classe de gestion du parametrage de SIFT
//---------------------------------------------------------
class ParametrageSift{
protected:
  int    _firstoctave;
  int    _nboctaves;
  int    _nblevels;
  /*float*/ double  _threshold;
  /*float*/ double  _edgeThreshold;
  int    _nodescr;
  int    _noorient;
  int    _savegss;
  int    _verbose;
  int    _binary;
  //int    _haveKeypoints;
  std::vector<T3D<float> > _KeyPointsADecrire;//cas ou l'on a une liste de points que l'on souhaite decrire
  int _siftdense; //pas pour faire du sift dense...
  /*float*/ double _siftdensescalefixe;
  /*signed char*/signed int  _signe;//signe < > : 0 > et < ; -1 < ; 1 >
  /*float*/ double _sigman;
  //  float _sigma0;
  bool _demidescripteur;

  
  




  //	      const float sigman,
  //	      const float sigma0,

  //	      VL::PgmBuffer & buffer,
  //	      const stringpp & outputFilename,
  //	      const stringpp & keypointsFilename,
  //	      const int xmin, const int ymin, const int xmax, const int ymax,
  //	      const double xcrop0, const double ycrop0,
  //	      std::ofstream & fkey,
  //	      unsigned int & nbkp

public:
  ParametrageSift(){defaut();}
  virtual ~ParametrageSift(){};

  //Parametrage par defaut
  void defaut(){
    _firstoctave    = -1 ;
    _nboctaves        = -1 ;
    _nblevels         = 3 ;
    _threshold      = 0.04f / _nblevels / 2.0f ;
    _edgeThreshold  = 10.0f;
    _nodescr        = 0 ;
    _noorient       = 0 ;
    _savegss        = 0 ;
    _binary         = 0 ;//plus vraiment d'utilité actuellement... ??
    //    _haveKeypoints  = 0 ;
    _signe          = 0 ;//signe < > : 0 > et < ; 1 < ; 2 >
    _sigman         = .5 ;
    //    _sigma0         = 1.6 * powf(2.0f, 1.0f / _nblevels) ;
    _siftdense=0;
    _siftdensescalefixe=-1.0;//si <0 echelle=premier niveau de la premiere octave // ancienne version 1.0;//echelle par defaut
    _demidescripteur=false;
    _verbose=false;
 
  }

  //Au cas ou l'on fait du sift dense : on ne travaille qu'a une seule echelle, et donc une seule octave
  //On peut donc revenir sur certains parametres...
  void update_cas_siftdense(){
    if(issiftdense()){return;}
    if(issiftdensescalefixe()){
      nboctaves()=1;
      firstoctave()=floor(log2((double(siftdensescalefixe()))/(double(sigma0())))-1.0/(double)nblevels());
    }
    nboctaves()=1;
  }
      
  //Parametres SIFT
  const int firstoctave() const{return _firstoctave;}
  const int nboctaves() const{return _nboctaves;}
  const int nblevels() const{return _nblevels;}
  const /*float*/double threshold() const{return _threshold;}
  const /*float*/ double edgeThreshold() const{return _edgeThreshold;}
  const int nodescr() const{return _nodescr;}
  const int noorient() const{return _noorient;}
  const int savegss() const{return _savegss;}
  const int binary() const{return _binary;}
  const int haveKeypoints() const{return /*_haveKeypoints*/_KeyPointsADecrire.size();}
  const /*signed char*/ signed int signe() const{return _signe;}//signe < > : 0 > et < ; 1 < ; 2 >
  const /*float*/ double sigman() const{return _sigman;}
  const /*float*/ double sigma0() const{return 1.6 * powf(2.0f, 1.0f / _nblevels);}
  const std::vector<T3D<float> > & KeyPointsADecrire() const{return _KeyPointsADecrire;}
  const int siftdense() const{return _siftdense;}
  const bool issiftdense() const{return _siftdense>0;}
  const bool issiftdensescalefixe() const{return _siftdensescalefixe>0;}
  const /*float*/ double siftdensescalefixe() const{return _siftdensescalefixe;}
  const bool demidescripteur() const{return _demidescripteur;}

  //Parametres SIFT
  int & firstoctave() {return _firstoctave;}
  int & nboctaves() {return _nboctaves;}
  int & nblevels() {return _nblevels;}
  /*float*/double & threshold() {return _threshold;}
  /*float*/ double & edgeThreshold() {return _edgeThreshold;}
  int & nodescr() {return _nodescr;}
  int & noorient() {return _noorient;}
  int & savegss() {return _savegss;}
  int & binary() {return _binary;}
  //  int & haveKeypoints() {return _haveKeypoints;}
  /*signed char*/signed int & signe() {return _signe;}//signe < > : 0 > et < ; -1 < ; 1 >
  std::vector<T3D<float> > & KeyPointsADecrire() {return _KeyPointsADecrire;}
  int & siftdense() {return _siftdense;}
  /*float*/ double & siftdensescalefixe() {return _siftdensescalefixe;}
  bool & demidescripteur() {return _demidescripteur;}

  const bool lireKeyPointsADecrire(const std::string & keypointsFilename);
  const bool lireKeyPointsADecrire_formatbin(const std::string & keypointsFilename);
  const bool lireKeyPointsADecrire_formattxt(const std::string & keypointsFilename);
  const bool lireKeyPointsADecrire_formatbin(std::ifstream & fkp);
  const bool lireKeyPointsADecrire_formattxt(std::ifstream & fkp);


};





//=========================================================
//ParametrageSiftDriver
//
//Classe de gestion des differents parametres utilises lors de la mise en oeuvre de sift
//  - ParametrageDallage : parametres lies au dallage des images
//  - ParametragePretraitement : parametres lies a divers pretraitements (choix du canal, gradient ? , reetalement)
//  - ParametrageSift : parametres lies directement a sift
//  - ParametrageSiftDriver : regroupement des parametres precedents
//---------------------------------------------------------
class ParametrageSiftDriver{
protected:
  ParametragePretraitement _param_pretraitement;
  ParametrageDallage _param_dallage;
  ParametrageSift _param_sift;
  int _verbose;
public:
  ParametrageSiftDriver(){defaut();}
  virtual ~ParametrageSiftDriver(){}

  void defaut(){
    _param_pretraitement.defaut();
    _param_dallage.defaut();
    _param_sift.defaut();
    _verbose        = 0 ;
  }
  
  const ParametragePretraitement & pretraitement() const{return _param_pretraitement;}
  const ParametrageDallage & dallage() const{return _param_dallage;}
  const ParametrageSift & sift() const{return _param_sift;}
  const int verbose() const{return _verbose;}

  ParametragePretraitement & pretraitement() {return _param_pretraitement;}
  ParametrageDallage & dallage() {return _param_dallage;}
  ParametrageSift & sift() {return _param_sift;}
  int & verbose() {return _verbose;}


  const bool LectureFichierParametres(std::string & adresse_fichier_parametres);
  const bool LectureFichierParametres(std::ifstream & ficparam);

};




#endif
