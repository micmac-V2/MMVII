#ifndef _KEYPOINT_HPP_
#define _KEYPOINT_HPP_

#include <vector>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>


#include "t2d.hpp"
#include "t3d.hpp"
#include "trectangle.hpp"
#ifdef _AVECMATIS_
//#warning avecmatis
//#include "matis/tpoint2d.hpp"
#include "matis/tchaine.hpp"
//#include "matis/trectangle.hpp"
#include "matis/ttimage.hpp"
#else
//#warning sansmatis
#endif

//#include "OutilsString.h"
#include "outilsstring/Stringpp.hpp"

//------------------------------------------------------------------------
//Classe de gestion des points d'intérêts
//------------------------------------------------------------------------
class KeyPoint{
public:
  typedef std::vector<KeyPoint> KeyPoints;
 
  ~KeyPoint(){};
  //Constructeurs
  KeyPoint(const int tailledesc=128);
  KeyPoint(const double x, const double y, const int tailledesc=128);
  KeyPoint(const double x, const double y, const double sigma, const double orientation, const int tailledesc=128);
  KeyPoint(const double x, const double y, const double sigma, const double orientation, const bool signe, const int tailledesc=128);
  KeyPoint(const double x, const double y, const double sigma, const double orientation, const std::vector<unsigned char> & desc);
  KeyPoint(const double x, const double y, const double sigma, const double orientation, const bool signe, const std::vector<unsigned char> & desc);
  KeyPoint(const T2D<double> & xy, const int tailledesc=128);
  KeyPoint(const T2D<double> & xy, const double sigma, const double orientation, const int tailledesc=128);
  KeyPoint(const T2D<double> & xy, const double sigma, const double orientation, const bool signe, const int tailledesc=128);
  KeyPoint(const T2D<double> & xy, const double sigma, const double orientation, const std::vector<unsigned char> & desc);
  KeyPoint(const T2D<double> & xy, const double sigma, const double orientation, const bool signe, const std::vector<unsigned char> & desc);

  //Fonctions d'acces aux donnees
  //  const T2D<double> & xy(){return T2D<double>(_x,_y);}//ne peut fonctionner que si l'on stocke les coordonnées du point sous la forme d'un objet T2D
  const int dim() const {return _desc.size();} 
  const int tailledesc() const{return _desc.size();}
  const T2D<double> xy() const {return T2D<double>(_x,_y);}
  const T2D<double> xy() {return T2D<double>(_x,_y);}
  void xy(const T2D<double> ptxy){_x=ptxy.X(),_y=ptxy.Y();}
  const KeyPoint sansdesc() const{return KeyPoint(_x,_y,_sigma,_orientation,_signe,0);}
  double & x(){return _x;} 
  const double x() const {return _x;} 
  double & y(){return _y;}
  const double y() const {return _y;}
  double & X(){return _x;} 
  const double X() const {return _x;} 
  double & Y(){return _y;}
  const double Y() const {return _y;}
  double & sigma(){return _sigma;} 
  const double sigma() const {return _sigma;} 
  double & orientation(){return _orientation;}
  const double orientation() const {return _orientation;}
  double & angle(){return _orientation;}
  const double angle() const {return _orientation;}
  const bool signe() const {return _signe;}
  bool & signe(){return _signe;}
  const unsigned int taillevoisinage() const{//cf calcul orientation siftpp...
    //    double winFactor=1.5; double sigmaw=winFactor*sigma ; unsigned int W=(unsigned int)ceil(3.0*sigmaw); return W;
    return (unsigned int)ceil(4.5*_sigma);
  }
  const std::vector<unsigned char> desc() const {return _desc;}
  std::vector<unsigned char> & desc() {return _desc;}
  unsigned char & operator[](const int i){return _desc[i];}
  const unsigned char operator[](const int i) const {return _desc[i];}  
  const bool operator == (const KeyPoint & kp);
  const bool operator != (const KeyPoint & kp);

  const KeyPoint & operator = (const KeyPoint & kp){
    _x=kp._x;
    _y=kp._y;
    _sigma=kp._sigma;
    _orientation=kp._orientation;
    _signe=kp._signe;
    _desc=kp._desc;
    return (*this);
  }
  
  inline static const double distance(const KeyPoint & kp1, const KeyPoint & kp2);
  inline static const double distanceeuclidienne(const KeyPoint & kp1, const KeyPoint & kp2);
  inline static const double distanceeuclidiennecarree(const KeyPoint & kp1, const KeyPoint & kp2);
  inline static const double distanceL1(const KeyPoint & kp1, const KeyPoint & kp2);
  
  //permet de normaliser le vecteur descripteur
  const std::vector<double> normalise(const double newnorme=1.0) const ;
  //calcule la norme du descripteur
  const double norme() const;
  //on ne prend plus en compte le sens du gradient mais seulement sa direction(approximation...)
  void DemiDescripteur(const int nbins=8, const bool changetailledesc=true);

  //Fonctions de lecture et d'écriture des points d'intérêt et de leur descripteur dans des fichiers .key générés par l'exécutable de Lowe
  void lireKeyLowe(std::ifstream & f, const int tailledesc);
  void lireKeyLowe(std::ifstream & f);
  void ecrireKeyLowe(std::ofstream & f) const;

//Fonctions de lecture et d'écriture des points d'intérêt et de leur descripteur dans des fichiers binaires
  void lireBin(std::ifstream & f, const int tailledesc);
  void lireBin(std::ifstream & f);
  void ecrireBin(std::ofstream & f) const;

  //Fonctions de lecture et d'écriture des points d'intérêt et de leur descripteur dans des fichiers .key générés par l'exécutable de l'INRIA
  void lireKeyINRIA(std::ifstream & f, const int tailledesc);
  void lireKeyINRIA(std::ifstream & f);
  void ecrireKeyINRIA(std::ofstream & f) const;

  //Fonctions de lecture et d'écriture des points d'intérêt et de leur descripteur dans des fichiers .key générés par SURF
  void lireKeySURF(std::ifstream & f, const int tailledesc);
  void lireKeySURF(std::ifstream & f);

  //Fonctions de lecture et d'écriture des points d'intérêt et de leur descripteur dans des fichiers .xml du type autopano-sift
  //  void lireKeyxml(std::ifstream & f);
  void ecrireKeyxml(std::ofstream & f) const;


  //Fonctions de lecture et d'écriture des points d'intérêt et de leur descripteur dans des fichiers .key générés par l'exécutable siftpp
  void lireKeysiftpp(std::ifstream & f);
  void ecrireKeysiftpp(std::ofstream & f) const;


  //Fonctions de lecture et d'écriture du point d'intérêt et de son descripteur dans des fichiers .loc et .desc
  void ecrireLoc(std::ofstream & f) const;
  void ecrireDesc(std::ofstream & f) const;

  void ecrireLocDesc(std::ofstream & floc, std::ofstream & fdesc) const;
  void lireLocDesc(std::ifstream & floc, std::ifstream & fdesc, int tailledesc);
  void lireLocDesc(std::ifstream & floc, std::ifstream & fdesc);

  typedef std::vector<unsigned char>::iterator iterator;
  typedef std::vector<unsigned char>::const_iterator const_iterator;
  iterator begin () {return _desc.begin();}
  const_iterator begin () const {return _desc.begin() ;}	 
  iterator end () {return _desc.end();}
  const_iterator end () const {return _desc.end();}  

protected:
  double _x,_y,_sigma,_orientation;
  bool _signe;
  std::vector<unsigned char> _desc;

  void initdesc(const int tailledesc=128){_desc=std::vector<unsigned char>(tailledesc);}//crée un vecteur descripteur "vide" de taille tailledesc
};



//------------------------------------------------------------------------
//Manipulation de fichier de points d'intérêts
//------------------------------------------------------------------------
namespace OutilsIOKeyPoint{
//Traduction d'un fichier sift de Lowe <->fichier binaire
  void KeyLowe2KeyBin(const stringpp & adresse_keylowe, const stringpp & adresse_keybin);
  void KeyBin2KeyLowe(const stringpp & adresse_keybin, const stringpp & adresse_keylowe);
  void KeyINRIA2KeyBin(const stringpp & adresse_keyinria, const stringpp & adresse_keybin);
  void KeyBin2KeyINRIA(const stringpp & adresse_keybin, const stringpp & adresse_keyinria);
  void KeyLowe2KeyINRIA(const stringpp & adresse_keylowe, const stringpp & adresse_keyinria);
  void KeyINRIA2KeyLowe(const stringpp & adresse_keyinria, const stringpp & adresse_keylowe);
  void Keykoen2KeyBin(const stringpp & adresse_keykoen, const stringpp & adresse_keybin);
  void KeyBin2Keykoen(const stringpp & adresse_keybin, const stringpp & adresse_keykoen, const bool avecdesc=true);
//Traduction binaire <-> xml
//  void KeyBin2Keyxml(const stringpp & adresse_keybin, const stringpp & adresse_keyxml);
//  void Keyxml2KeyBin(const stringpp & adresse_keyxml, const stringpp & adresse_keybin);
//Traduction d'un fichier siftpp en un fichier sift de Lowe
  void Siftpp2KeyLowe(const stringpp & adresse_siftpp, const stringpp & adresse_keylowe);
// Traduction d'un fichier .key (généré par l'exécutable de Lowe) en un fichier .loc et en un fichier .desc (la version _card renvoie en plus le nombre de points...)
  void KeyLowe2LocDesc(const stringpp & adresse_key, const stringpp & adresse_loc, const stringpp & adresse_desc);
  void KeyLowe2LocDesc(const stringpp & adresse);
  void KeyLowe2LocDesc_card(const stringpp & adresse_key, const stringpp & adresse_loc, const stringpp & adresse_desc, int & nbpts);
  void KeyLowe2LocDesc_card(const stringpp & adresse, int & nbpts);
  void KeyLowe2LocDesc_card(const stringpp & adresse_key, const stringpp & adresse_loc, const stringpp & adresse_desc, int & nbpts, int & tailledesc);
  void KeyLowe2LocDesc_card(const stringpp & adresse, int & nbpts, int & tailledesc);

// Connaitre le nombre de points contenus dans le fichier
// Connaitre la taille du descripteur
  const int card_KeyLowe(const stringpp & adresse_key);
  const int tailledesc_KeyLowe(const stringpp & adresse_key);
  void entete_KeyLowe(const stringpp & adresse_key, int & nbpts, int & tailledesc);
  const int card_KeyBin(const stringpp & adresse_key);
  const int card_KeyBin(std::ifstream & f_key);
  const int tailledesc_KeyBin(const stringpp & adresse_key);
  const int tailledesc_KeyBin(std::ifstream  & f_key);
  void entete_KeyBin(const stringpp & adresse_key, int & nbpts, int & tailledesc);
  void entete_KeyBin(std::ifstream & f_key, int & nbpts, int & tailledesc);
  const int card_KeyINRIA(const stringpp & adresse_key);
  const int tailledesc_KeyINRIA(const stringpp & adresse_key);
  void entete_KeyINRIA(const stringpp & adresse_key, int & nbpts, int & tailledesc);



  const int card_Loc(const stringpp & adresse_loc); 
  const int card_Desc(const stringpp & adresse_pts, const int tailledesc=128); 
  const int card_LocDesc(const stringpp & adresse_loc, const stringpp & adresse_desc, int tailledesc=128); 
  const int card_LocDesc(const stringpp & adresse, const int tailledesc=128); 

//Traduction d'un fichier .loc et d'un fichier .desc en un fichier .key (du type de ceux générés par l'exécutable de Lowe)(la version _card renvoie en plus le nombre de points)
  void LocDesc2KeyLowe(const stringpp & adresse_loc, const stringpp & adresse_desc, const stringpp & adresse_key, int tailledesc=128);
  void LocDesc2KeyLowe(const stringpp & adresse, const int tailledesc=128);
  void LocDesc2KeyLowe_card(const stringpp & adresse_loc, const stringpp & adresse_desc, const stringpp & adresse_key, int & nbpts, const int tailledesc=128);
  void LocDesc2KeyLowe_card(const stringpp & adresse, int & nbpts, const int tailledesc=128);

//Chargement d'un fichier .key (généré par l'exécutable de Lowe). On obtient en sortie un vector de points d'intérêt (objets KeyPoint)
  const bool ChargerKeyLowe(const stringpp & adresse_key, std::vector<KeyPoint> & kp);
  const bool ChargerKeyLowe(const stringpp & adresse_key, std::vector<T2D<double> > & kp);
  const bool ChargerLoc(const stringpp & adresse_loc, std::vector<T2D<double> > & kp);

  const bool ChargerKeyBin(const stringpp & adresse_key, std::vector<KeyPoint> & kp);
  const bool ChargerKeyBin(const stringpp & adresse_key, std::vector<T2D<double> > & kp);
  const bool ChargerKeyBin(const stringpp & adresse_key, std::vector<KeyPoint> & kp, const double sigmamin, const double sigmamax);

  const bool ChargerKeyINRIA(const stringpp & adresse_key, std::vector<KeyPoint> & kp);
  const bool ChargerKeyINRIA(const stringpp & adresse_key, std::vector<T2D<double> > & kp);

  const bool ChargerKeyxml(const stringpp & adresse_key, std::vector<KeyPoint> & kp);
  const bool ChargerKeyxml(const stringpp & adresse_key, std::vector<T2D<double> > & kp);

//Chargement d'un fichier .loc et d'un fichier .desc . On obtient en sortie un vector de points d'intérêt (objets KeyPoint)
  const bool ChargerLocDesc(const stringpp & adresse_loc, const stringpp & adresse_desc, std::vector<KeyPoint> & kp, const int tailledesc=128);
  const bool ChargerLocDesc(const stringpp & adresse, std::vector<KeyPoint> & kp, const int tailledesc=128);

//Ecriture d'une liste de points d'intérêt et leur descripteur SIFT (objets KeyPoint)
  const bool EcrireKeyLowe(const stringpp & adresse_key, const std::vector<KeyPoint> & kp);
  const bool EcrireKeyINRIA(const stringpp & adresse_key, const std::vector<KeyPoint> & kp);
  const bool EcrireKeyxml(const stringpp & adresse_key, const std::vector<KeyPoint> & kp);
  const bool EcrireKeyBin(const stringpp & adresse_key, const std::vector<KeyPoint> & kp);

//Ecriture d'une liste de points d'intérêt (uniquement leurs coordonnées...) donnée sous la forme d'un vector d'objets KeyPoint dans un fichier .loc
  const bool EcrireLoc(const stringpp & adresse_loc, const std::vector<KeyPoint> & kp);

//Ecriture d'une liste de points d'intérêt (uniquement leur descripteur SIFT...) donnée sous la forme d'un vector d'objets KeyPoint dans un fichier .desc
  const bool EcrireDesc(const stringpp & adresse_desc, const std::vector<KeyPoint> & kp);

//Ecriture d'une liste de points d'intérêt donnés sous la forme d'un vector d'objets KeyPoint dans un fichier .loc (pour leurs coordonnées) et dans un fichier .desc (pour leur descripteur SIFT)
  const bool EcrireLocDesc(const stringpp & adresse_loc, const stringpp & adresse_desc, const std::vector<KeyPoint> & kp);
  const bool EcrireLocDesc(const stringpp & adresse, const std::vector<KeyPoint> & kp);
};

namespace OutilsKeyPoint{
//Calcul de la distance entre les descripteurs SIFT de deux points d'intérêt
  const double distance(const KeyPoint & kp1, const KeyPoint & kp2);
  const double distanceeuclidienne(const KeyPoint & kp1, const KeyPoint & kp2);
  const double distanceeuclidiennecarree(const KeyPoint & kp1, const KeyPoint & kp2);
  const double distanceL1(const KeyPoint & kp1, const KeyPoint & kp2);


//Sélection des points d'intérêt compris dans une fenêtre
  void selection_fenetre(const std::vector<KeyPoint> & kp, const TRectangle2D<float> & fenetre, std::vector<KeyPoint> & kpsel, bool strict=false);
  void selection_fenetre(std::vector<KeyPoint> & kp, const TRectangle2D<float> & fenetre, bool strict=false);

//Descripteur moyen de tous les descripteurs
  void DescripteurMoyen(const std::vector<KeyPoint> & vkp, KeyPoint & kpm);
  void DescripteurMoyenGrille(const std::vector<KeyPoint> & vkp, int nc, int nl, double dc, double dl, std::vector<KeyPoint> & vkpm);

#ifdef _AVECMATIS_
//Filtrage des points hors du masque
  void FiltrageMasque(const std::vector<KeyPoint> & kp, const TTImage<bool> & masque, std::vector<KeyPoint> & kpsel);
  void FiltrageMasque(std::vector<KeyPoint> & kp, const TTImage<bool> & masque);
#endif

//typedef
  typedef std::vector<KeyPoint> KeyPoints;
  


  

  

};



namespace OutilsIOKeyPoint{
  //====================================================
  //Classes Flux de points d'interet
  //----------------------------------------------------

  //Classe de base
  class FluxSortieKeypoint{
  public:
    FluxSortieKeypoint(){}
    virtual ~FluxSortieKeypoint(){};
    virtual const bool operator << (const KeyPoint & kp) = 0;
    virtual const unsigned int size() const = 0; 
  };
  

  //Flux dans un fichier binaire
  class FluxSortieKeypoint_fichierbinaire : FluxSortieKeypoint{
  protected:
    std::ofstream & _out;
    unsigned int _nbkp;
    unsigned int _dim;

    void commencer(){
      //On prépare les champs d'entete, qui seront rempli a la fin...
      _out.write((const char*)&_nbkp,sizeof(unsigned int));
      _out.write((const char*)&_dim,sizeof(unsigned int));
    }

  public:
    FluxSortieKeypoint_fichierbinaire(std::ofstream & fout):_out(fout),_nbkp(0),_dim(128){commencer();}
    FluxSortieKeypoint_fichierbinaire(std::ofstream & fout, const unsigned int dim):_out(fout),_nbkp(0),_dim(dim){commencer();}
    virtual ~FluxSortieKeypoint_fichierbinaire(){}
    virtual const bool operator << (const KeyPoint & kp){
      kp.ecrireBin(_out);
      _dim=kp.dim();
      _nbkp++;
	  return true;
    }
    void finir(){
      //On peut maintenant ecrire le nombre de points dans le fichier
      _out.seekp(0,std::ios::beg);
      _out.write((const char*)&_nbkp,sizeof(unsigned int));
      _out.write((const char*)&_dim,sizeof(unsigned int));
      _out.close();
    }
    void close(){finir();}
    virtual const unsigned int size() const{return _nbkp;}
    virtual const unsigned int dim() const{return _dim;}
  };


  //Flux dans un vector
  class FluxSortieKeypoint_vector : FluxSortieKeypoint{
  protected :
    std::vector<KeyPoint> & _out;
  public:
  FluxSortieKeypoint_vector(std::vector<KeyPoint> & out):_out(out){}
    virtual ~FluxSortieKeypoint_vector(){}
    virtual const bool operator << (const KeyPoint & kp){
      _out.push_back(kp);
      return true;
    }  
    virtual const unsigned int size() const{return _out.size();}
    virtual const unsigned int dim() const{return (_out).begin()->dim();}
  };


  //====================================================
  //Classes Flux de lecture de points d'interet
  //----------------------------------------------------
  //Classe de base
  class FluxLectureKeypoint{
  public:
    FluxLectureKeypoint(){}
    virtual ~FluxLectureKeypoint(){};
    virtual const bool operator >> (KeyPoint & kp) = 0;
    virtual const unsigned int size() const = 0; 
    virtual const unsigned int dim() const = 0; 
  };

  //Flux de lecture dans un fichier binaire
  class FluxLectureKeypoint_fichierbinaire : FluxLectureKeypoint{
  protected:
    std::ifstream & _in;
    unsigned int _nbkp;
    unsigned int _dim;
    unsigned int _pos;

    void commencer(){
      _in.read((char*)&_nbkp,sizeof(unsigned int));
      _in.read((char*)&_dim,sizeof(unsigned int));
      _pos=0;
    }
    void commencer(const int pos){
      _in.read((char*)&_nbkp,sizeof(unsigned int));
      _in.read((char*)&_dim,sizeof(unsigned int));
      atteindre(pos);
    }

    void atteindre(const unsigned int pos){
      if(_pos<pos){
	_in.seekg(0,std::ios::beg);
	unsigned int tmp;
	_in.read((char*)&tmp,sizeof(unsigned int));
	_in.read((char*)&tmp,sizeof(unsigned int));
	KeyPoint kptmp(_dim);
	_pos=0;
	for(unsigned int n=0;n<pos;n++){operator>>(kptmp);}
      }else{
	KeyPoint kptmp(_dim);
	for(unsigned int n=0;n<(pos-_pos);n++){operator>>(kptmp);}
      }
    }

  public:
    FluxLectureKeypoint_fichierbinaire(std::ifstream & fin):_in(fin),_nbkp(0),_dim(128),_pos(0){commencer();}
  FluxLectureKeypoint_fichierbinaire(std::ifstream & fin, const int pos):_in(fin),_nbkp(0),_dim(128),_pos(pos){commencer(pos);}
    virtual ~FluxLectureKeypoint_fichierbinaire(){}
    virtual const bool operator >> (KeyPoint & kp){
      if(_pos>=size()){
	return false;
      }
      kp.lireBin(_in);
      return true;
    }
    void finir(){}
    void close(){}
    virtual const unsigned int size() const{return _nbkp;}
    virtual const unsigned int dim() const{return _dim;}
  };


  //Flux de lecture dans un fichier texte KeyLowe
  class FluxLectureKeypoint_fichierKeyLowe : FluxLectureKeypoint{
  protected:
    std::ifstream & _in;
    unsigned int _nbkp;
    unsigned int _dim;
    unsigned int _pos;

    void commencer(){
      _in>>_nbkp;
      _in>>_dim;
      _pos=0;
    }
    void commencer(const int pos){
      _in>>_nbkp;
      _in>>_dim;
      atteindre(pos);
    }

    void atteindre(const unsigned int pos){
      if(_pos<pos){
	_in.seekg(0,std::ios::beg);
	unsigned int tmp;
	_in>>tmp>>tmp;
	KeyPoint kptmp(_dim);
	_pos=0;
	for(unsigned int n=0;n<pos;n++){operator>>(kptmp);}
      }else{
	KeyPoint kptmp(_dim);
	for(unsigned int n=0;n<(pos-_pos);n++){operator>>(kptmp);}
      }
    }

  public:
    FluxLectureKeypoint_fichierKeyLowe(std::ifstream & fin):_in(fin),_nbkp(0),_dim(128),_pos(0){commencer();}
    FluxLectureKeypoint_fichierKeyLowe(std::ifstream & fin, const int pos):_in(fin),_nbkp(0),_dim(128),_pos(pos){commencer(pos);}
    virtual ~FluxLectureKeypoint_fichierKeyLowe(){}
    virtual const bool operator >> (KeyPoint & kp){
      if(_pos>=size()){
	return false;
      }
      kp.lireKeyLowe(_in);
      return true;
    }
    void finir(){}
    void close(){}
    virtual const unsigned int size() const{return _nbkp;}
    virtual const unsigned int dim() const{return _dim;}
  };


  //Flux de lecture dans un vector
  class FluxLectureKeypoint_vector : FluxLectureKeypoint{
  protected:
    const std::vector<KeyPoint> & _in;
    unsigned int _dim;
    unsigned int _pos;

    void commencer(){
      _dim=_in.front().dim();
      _pos=0;
    }
    void commencer(const int pos){
      _dim=_in.front().dim();
      atteindre(pos);
    }

    void atteindre(const int pos){
      _pos=pos;
    }

  public:
    FluxLectureKeypoint_vector(const std::vector<KeyPoint> & vin):_in(vin),_pos(0){commencer();}
    FluxLectureKeypoint_vector(const std::vector<KeyPoint> & vin, const int pos):_in(vin),_pos(pos){commencer(pos);}
    virtual ~FluxLectureKeypoint_vector(){}
    virtual const bool operator >> (KeyPoint & kp){
      if(_pos>=size()){
	return false;
      }
      kp=_in[_pos];
      return true;
    }
    void finir(){}
    void close(){}
    virtual const unsigned int size() const{return _in.size();}
    virtual const unsigned int dim() const{return _dim;}
  };



};



#endif
