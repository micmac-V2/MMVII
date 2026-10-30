#ifndef RESULTMEC_H
#define RESULTMEC_H

#include <vector>
#include <map>
#include <fstream>

#include "t2d.hpp"
#include "trectangle.hpp"

#ifdef _AVECMATIS_
//#include "matis/tchaine.hpp"
//#include "matis/tpoint2d.hpp"
//#include "matis/trectangle.hpp"
//#include "matis/macro_sortiemessage.hpp"
#include "matis/orientation.hpp"
#include "matis/ttimage.hpp"
#else
//#include "tpoint2d.hpp"
//#include "trectangle.hpp"
#define SORTIEMESSAGE(a) std::cout<<a;
#endif

#include "outilsstring/Stringpp.hpp"
#include <vector>


//using namespace std;


//#define ES_RESULT_BINAIRE

//------------------------------------------------------------------
//Classe de manipulation des resultats de mise en correspondance : points homologues entre 2 images
//------------------------------------------------------------------
template<class T>
class ResultMEC{
 public:
  //Constructeurs
  ResultMEC():_nom1(""),_nom2(""){};
  ResultMEC(const stringpp & adresseresult){DecombineNom(adresseresult,_nom1,_nom2);Lire(adresseresult);}
  ResultMEC(const stringpp & adresseresult, const stringpp & nom1, const stringpp & nom2):_nom1(nom1),_nom2(nom2){Lire(adresseresult);}
  ResultMEC(const std::vector<T> & pts1, const std::vector<T> & pts2):_nom1(""),_nom2(""),_pts1(pts1),_pts2(pts2){};
  ResultMEC(const stringpp & nom1, const stringpp & nom2):_nom1(nom1),_nom2(nom2){};
  ResultMEC(const stringpp & nom1, const std::vector<T> & pts1, const stringpp & nom2, const std::vector<T> & pts2):_nom1(nom1),_nom2(nom2),_pts1(pts1),_pts2(pts2){};
  ResultMEC(const unsigned int nbelts):_nom1(""),_nom2(""),_pts1(nbelts),_pts2(nbelts){};
  ~ResultMEC(){};

  //typedef
  typedef T TypeObjet;

  //swap
  void swap(){std::swap(_pts1,_pts2);}
  
  //Outils pour manipuler les noms de fichier image1 image2 <-> image1.-.image2
  void DecombineNom(const stringpp & adresseresult, stringpp & nom1, stringpp & nom2) const{stringpp::DeCombineNom(adresseresult,nom1,nom2);};
  void DecombineNom(const stringpp & adresseresult){DecombineNom(adresseresult,_nom1,_nom2);}
  const stringpp CombineNom(const stringpp & repertoire, const stringpp & nom1, const stringpp & nom2) const{return stringpp::AdressefromDirFile(repertoire,CombineNom(nom1,nom2));}
  const stringpp CombineNom(const stringpp & nom1, const stringpp & nom2) const{return stringpp::CombineNom(nom1,nom2,stringpp("result"));}
  const stringpp CombineNom() const{return CombineNom(_nom1,_nom2);}
  const stringpp CombineNom(const stringpp & repertoire) const{return CombineNom(repertoire,_nom1,_nom2);}
  
  //Renvoie le nom des images (si elles en ont)
  const stringpp nom1() const{return _nom1;}
  stringpp & nom1(){return _nom1;}
  const stringpp nom2() const{return _nom2;}
  stringpp & nom2(){return _nom2;}
  
  //Fonctions de Lecture/Ecriture avec des noms de fichier du type image1.-.image2.result
  const bool Charger(const stringpp & adresseresult){clear();DecombineNom(adresseresult,_nom1,_nom2);return Lire(adresseresult);}
  const bool ChargerTexte(const stringpp & adresseresult){clear();DecombineNom(adresseresult,_nom1,_nom2);return LireTexte(adresseresult);}
  const bool ChargerBinaire(const stringpp & adresseresult){clear();DecombineNom(adresseresult,_nom1,_nom2);return LireBinaire(adresseresult);}
  void Sauver(const stringpp & adresserepertoireresult) const{Ecrire(CombineNom(adresserepertoireresult));}
  void SauverTexte(const stringpp & adresserepertoireresult) const{EcrireTexte(CombineNom(adresserepertoireresult));}
  void SauverTexte(const stringpp & adresserepertoireresult, const int prec) const{EcrireTexte(CombineNom(adresserepertoireresult),prec);}
  void SauverBinaire(const stringpp & adresserepertoireresult) const{EcrireBinaire(CombineNom(adresserepertoireresult));}

  //Fonctions de simple Lecture/Ecriture (on ne s'interesse qu'au contenu du fichier)
  const bool Lire(const stringpp & adresseresult);
  const bool LireBinaire(const stringpp & adresseresult);
  const bool LireTexte(const stringpp & adresseresult);
  const bool LireXML(const stringpp & adresseresult);
  const bool Lire(std::ifstream & fresult);
  const bool LireBinaire(std::ifstream & fresult);
  const bool LireTexte(std::ifstream & fresult);
  const bool Ecrire(const stringpp & adresseresult) const;
  const bool EcrireBinaire(const stringpp & adresseresult) const;
  const bool EcrireTexte(const stringpp & adresseresult) const;
  const bool EcrireTexte(const stringpp & adresseresult, const int prec) const;
  const bool Ecrire(std::ofstream & fresult) const;
  const bool EcrireBinaire(std::ofstream & fresult) const;
  const bool EcrireTexte(std::ofstream & fresult) const;
  const bool EcrireTexte(std::ofstream & fresult, const int prec) const;

  template<class U> static void Ecrire(std::ofstream & fresult, const T2D<U> & pt1, const T2D<U> & pt2);
  template<class U> static void EcrireBinaire(std::ofstream & fresult, const T2D<U> & pt1, const T2D<U> & pt2);
  template<class U> static void EcrireTexte(std::ofstream & fresult, const T2D<U> & pt1, const T2D<U> & pt2);
  template<class U> static void EcrireTexte(std::ofstream & fresult, const T2D<U> & pt1, const T2D<U> & pt2, const int prec);
  template<class U> static void Lire(std::ifstream & fresult, T2D<U> & pt1, T2D<U> & pt2);
  template<class U> static void LireBinaire(std::ifstream & fresult, T2D<U> & pt1, T2D<U> & pt2);
  template<class U> static void LireTexte(std::ifstream & fresult, T2D<U> & pt1, T2D<U> & pt2);
  
  
  //Fonctions de Lecture/Ecriture de pk0 (qui donnent des noms appropries aux fichiers)
  const bool ChargerPk0(const stringpp & adresse1, const stringpp & adresse2){
    clear();
    _nom1=adresse1.File();_nom1.SupprExtension();
    _nom2=adresse2.File();_nom2.SupprExtension();
    return LirePk0(adresse1,adresse2);
  }
  void SauverPk0(const stringpp & adresserepertoireresult, const int imin=0) const{
    EcrirePk0((stringpp::AdressefromDirFile(adresserepertoireresult,_nom1)).newSetExtension(".pk0"),(stringpp::AdressefromDirFile(adresserepertoireresult,_nom2)).newSetExtension(".pk0"),imin);
  }
  //Fonctions de Lecture/Ecriture de pk0 (on ne s'interesse qu'au contenu du fichier)
  const bool LirePk0(const stringpp & adresse1, const stringpp & adresse2);
  const bool LirePk0(std::ifstream & f1, std::ifstream & f2);
  const bool EcrirePk0(const stringpp & adresse1, const stringpp & adresse2, const int imin=0) const;
  const bool EcrirePK0(std::ofstream & f1, std::ofstream & f2, const int imin=0) const;
  

  const T pt1(const int id) const{return _pts1[id];}
  T & pt1(const int id){return _pts1[id];}
  const T pt2(const int id) const{return _pts2[id];}
  T & pt2(const int id){return _pts2[id];}
  const T2D<T> operator[] (const int id) const{return T2D<T>(_pts1[id],_pts2[id]);}
  void add(const T & pt1, const T & pt2){push_back(pt1,pt2);}
  void add(const T2D<T> & ptmec){push_back(ptmec);}
  void push_back(const T & pt1, const T & pt2){_pts1.push_back(pt1);_pts2.push_back(pt2);}
  void push_back(const T2D<T> & ptmec){_pts1.push_back(ptmec.x());_pts2.push_back(ptmec.y());}
  const std::vector<T> pts1() const{return _pts1;}
  std::vector<T> & pts1(){return _pts1;}
  const std::vector<T> pts2() const{return _pts2;}
  std::vector<T> & pts2(){return _pts2;}
  const unsigned int card() const{return _pts1.size();}//taille
  const unsigned int size() const{return _pts1.size();}//taille
  const bool empty() const{return _pts1.empty();}//vide ?
  const bool coherence1_2() const{return _pts1.size()==_pts2.size();}
  void resize(const int n){_pts1.resize(n); _pts2.resize(n);}

  const TRectangle2D<T> Emprise1() const;//Renvoie l'emprise des points dans la premiere image
  const TRectangle2D<T> Emprise2() const;//Renvoie l'emprise des points dans la seconde image
  
  //iterateurs pour aller plus vite
  typedef typename std::vector<T>::iterator iterator;
  typedef typename std::vector<T>::const_iterator const_iterator;
  iterator begin1(){return _pts1.begin();}
  iterator begin2(){return _pts2.begin();}
  const_iterator begin1() const {return _pts1.begin();}
  const_iterator begin2() const{return _pts2.begin();}
  iterator end1(){return _pts1.end();}
  iterator end2(){return _pts2.end();}
  const_iterator end1() const {return _pts1.end();}
  const_iterator end2() const{return _pts2.end();}
  
  //Nettoyage
  void clear(){_nom1=stringpp("");_nom2=stringpp("");_pts1.clear();_pts2.clear();}
#ifdef _AVECMATIS_
  //Pour conserver les mises en correspondances avec un masque...
  void FiltreMasque1(const TTImage<bool> & masque);
  void FiltreMasque2(const TTImage<bool> & masque);
  void FiltreMasque(const TTImage<bool> & masque1, const TTImage<bool> & masque2);
#endif

 protected:
  stringpp _nom1,_nom2;//Nom des images 1 et 2 (si elles en ont)
  std::vector<T> _pts1,_pts2;//Coordonnees des points homologues dans les images 1 et 2
};


//-------------------------------------------------------
//Fonctions de simple lecture/ecriture (On ne s'interesse qu'au contenu du fichier.)
//-------------------------------------------------------
template<class T> const bool ResultMEC<T>::Lire(const stringpp & adresseresult){
  std::ifstream fresult(adresseresult.c_str(),std::ios::in);
  if(!Lire(fresult)){
    SORTIEMESSAGE("Probleme en lecture du fichier "<<adresseresult<<"\n");
    return false;
  }
  fresult.close();
  return true;
}

template<class T> const bool ResultMEC<T>::LireBinaire(const stringpp & adresseresult){
  std::ifstream fresult(adresseresult.c_str(),std::ios::in | std::ios::binary);
  if(!LireBinaire(fresult)){
    SORTIEMESSAGE("Probleme en lecture du fichier "<<adresseresult<<"\n");
    return false;
  }
  fresult.close();
  return true;
}

template<class T> const bool ResultMEC<T>::LireTexte(const stringpp & adresseresult){
  std::ifstream fresult(adresseresult.c_str(),std::ios::in);
  if(!LireTexte(fresult)){
    SORTIEMESSAGE("Probleme en lecture du fichier "<<adresseresult<<"\n");
    return false;
  }
  fresult.close();
  return true;
}

template<class T> const bool ResultMEC<T>::LireXML(const stringpp & adresseresult){
  std::ifstream fxml(adresseresult.c_str(),std::ios::in);
  if(!fxml.good()){
    SORTIEMESSAGE("Fichier "<<adresseresult<<" introuvable\n");
    return false;
  }
  stringpp tmp;
  std::string::size_type pos=0;  
  float x1,y1,x2,y2;
  fxml>>tmp;
  while(true){
    fxml>>tmp;    
    //x1
    fxml>>tmp;
    pos=tmp.find_last_of("<");
    if(pos!=tmp.npos){tmp.erase(pos,tmp.length()-1);}
    pos=tmp.find_last_of(">");
    if(pos!=tmp.npos){tmp.erase(0,pos+1);}
    x1=atof(tmp.c_str());
    //y1
    fxml>>tmp;
    pos=tmp.find_last_of("<");
    if(pos!=tmp.npos){tmp.erase(pos,tmp.length()-1);}
    pos=tmp.find_last_of(">");
    if(pos!=tmp.npos){tmp.erase(0,pos+1);}
    y1=atof(tmp.c_str());
    //x2
    fxml>>tmp;
    pos=tmp.find_last_of("<");
    if(pos!=tmp.npos){tmp.erase(pos,tmp.length()-1);}
    pos=tmp.find_last_of(">");
    if(pos!=tmp.npos){tmp.erase(0,pos+1);}
    x2=atof(tmp.c_str());
    //y2
    fxml>>tmp;
    pos=tmp.find_last_of("<");
    if(pos!=tmp.npos){tmp.erase(pos,tmp.length()-1);}
    pos=tmp.find_last_of(">");
    if(pos!=tmp.npos){tmp.erase(0,pos+1);}
    y2=atof(tmp.c_str());
    fxml>>tmp;
    fxml>>tmp;
    if(fxml.eof()){break;}
    _pts1.push_back(T(x1,y1));
    _pts2.push_back(T(x2,y2));
  }
  return true;
}

template<class T> const bool ResultMEC<T>::Lire(std::ifstream & fresult){
#ifndef ESRESULTBINAIRE
  return LireTexte(fresult);
#else
  return LireBinaire(fresult);
#endif
}

template<class T> const bool ResultMEC<T>::LireBinaire(std::ifstream & fresult){
  //Passons a la version de micmac...
  if(!fresult.good())return false;
  int nbpts;
  int aa;
  fresult.read((char*)&(aa),sizeof(int));
  //std::cout<<aa<<std::endl;
  fresult.read((char*)&(nbpts),sizeof(int));
  //std::cout<<nbpts<<std::endl;
  fresult.read((char*)&(aa),sizeof(int));
  //std::cout<<aa<<std::endl;
  double a;
  T ptmp1,ptmp2;
  for(int n=0;n<nbpts;n++){
    fresult.read((char*)&(a),sizeof(double));
    //    std::cout<<"double avant"<<a<<std::endl;
    double dtmp;
    fresult.read((char*)&(dtmp),sizeof(double)); ptmp1.x()=dtmp;
    fresult.read((char*)&(dtmp),sizeof(double)); ptmp1.y()=dtmp;
    fresult.read((char*)&(dtmp),sizeof(double)); ptmp2.x()=dtmp;
    fresult.read((char*)&(dtmp),sizeof(double)); ptmp2.y()=dtmp;
    if(fresult.eof()){break;}
    _pts1.push_back(ptmp1);
    _pts2.push_back(ptmp2);
    fresult.read((char*)&(aa),sizeof(int));
    //    std::cout<<"int apres"<<aa<<std::endl;
  }
  
  return true;
  /*  
  if(!fresult.good())return false;
  T ptmp1,ptmp2;
  while(true){
    fresult.read((char*)&(ptmp1.x()),sizeof(float));
    fresult.read((char*)&(ptmp1.y()),sizeof(float));
    fresult.read((char*)&(ptmp2.x()),sizeof(float));
    fresult.read((char*)&(ptmp2.y()),sizeof(float));
    if(fresult.eof()){break;}
    _pts1.push_back(ptmp1);
    _pts2.push_back(ptmp2);
  }
  return true;
  */
}

template<class T> const bool ResultMEC<T>::LireTexte(std::ifstream & fresult){
  if(!fresult.good())return false;
  T ptmp1,ptmp2;
  while(true){
    fresult>>ptmp1.x();
    if(fresult.eof()){break;}//gestion du cas nouvelle ligne a la fin du fichier 
    fresult>>ptmp1.y()>>ptmp2.x()>>ptmp2.y();
    _pts1.push_back(ptmp1);
    _pts2.push_back(ptmp2);
    if(fresult.eof()){break;}
  }
  return true;
}

template<class T> const bool ResultMEC<T>::Ecrire(const stringpp & adresseresult) const{
  std::ofstream fresult(adresseresult.c_str(),std::ios::out);
  if(!fresult.good()){std::cout<<"Impossible d'ecrire dans "<<adresseresult<<std::endl; return false;}
  if(!Ecrire(fresult)){return false;}
  fresult.close();
  return true;
}

template<class T> const bool ResultMEC<T>::EcrireBinaire(const stringpp & adresseresult) const{
  std::ofstream fresult(adresseresult.c_str(),std::ios::out | std::ios::binary);
  if(!fresult.good()){std::cout<<"Impossible d'ecrire dans "<<adresseresult<<std::endl; return false;}
  if(!EcrireBinaire(fresult)){return false;}
  fresult.close();
  return true;
}

template<class T> const bool ResultMEC<T>::EcrireTexte(const stringpp & adresseresult) const{
  std::ofstream fresult(adresseresult.c_str(),std::ios::out);
  if(!fresult.good()){std::cout<<"Impossible d'ecrire dans "<<adresseresult<<std::endl; return false;}
  if(!EcrireTexte(fresult)){return false;}
  fresult.close();
  return true;
}

template<class T> const bool ResultMEC<T>::EcrireTexte(const stringpp & adresseresult, const int prec) const{
  if(prec<=0){return ResultMEC<T>::EcrireTexte(adresseresult);}
  std::ofstream fresult(adresseresult.c_str(),std::ios::out);
  if(!fresult.good()){std::cout<<"Impossible d'ecrire dans "<<adresseresult<<std::endl; return false;}
  if(!EcrireTexte(fresult,prec)){return false;}
  fresult.close();
  return true;
}

template<class T> const bool ResultMEC<T>::Ecrire(std::ofstream & fresult) const{
#ifndef ES_RESULT_BINAIRE
  return EcrireTexte(fresult);
#else
  return EcrireBinaire(fresult);
#endif
}

template<class T> const bool ResultMEC<T>::EcrireBinaire(std::ofstream & fresult) const{
  if(_pts1.size()!=_pts2.size()){return false;}
  if(_pts1.empty()){return true;}
  
  //Passons a la version de micmac...
  if(!fresult.good())return false;
  int nbpts=_pts1.size();
  int aa=2;
  fresult.write((const char*)&(aa),sizeof(int));
  fresult.write((const char*)&(nbpts),sizeof(int));
  fresult.write((const char*)&(aa),sizeof(int));
  double a=1;
  T ptmp1,ptmp2;
  for(int n=0;n<nbpts;n++){
    fresult.write((const char*)&(a),sizeof(double));
    double dtmp;
    dtmp=_pts1[n].X(); fresult.write((const char*)&(dtmp),sizeof(double)); 
    dtmp=_pts1[n].Y(); fresult.write((const char*)&(dtmp),sizeof(double)); 
    dtmp=_pts2[n].X(); fresult.write((const char*)&(dtmp),sizeof(double)); 
    dtmp=_pts2[n].Y(); fresult.write((const char*)&(dtmp),sizeof(double)); 
    fresult.write((const char*)&(aa),sizeof(int));
  }

  return true;


/*for(unsigned int n=0;n<_pts1.size();n++){
    fresult.write((const char*)&(_pts1[n].X()),sizeof(float));
    fresult.write((const char*)&(_pts1[n].Y()),sizeof(float));
    fresult.write((const char*)&(_pts2[n].X()),sizeof(float));
    fresult.write((const char*)&(_pts2[n].Y()),sizeof(float));
    }*/
/*
  const_iterator i1=begin1();
  const_iterator i2=begin2();
  for(;i1!=end1();++i1,++i2){
    fresult.write((const char*)&(*i1.X()),sizeof(float));
    fresult.write((const char*)&(*i1.Y()),sizeof(float));
    fresult.write((const char*)&(*i2.X()),sizeof(float));
    fresult.write((const char*)&(*i2.Y()),sizeof(float));
  }
  return true;
*/
}

template<class T> const bool ResultMEC<T>::EcrireTexte(std::ofstream & fresult, const int prec) const{
  if(_pts1.size()!=_pts2.size()){return false;}
  if(_pts1.empty()){return true;}
  if(prec<=0){return ResultMEC<T>::EcrireTexte(fresult);}
  /*  for(unsigned int n=0;n<_pts1.size();n++){
    fresult<<_pts1[n].X()<<"\t"<<_pts1[n].Y()<<"\t"<<_pts2[n].X()<<"\t"<<_pts2[n].Y()<<std::endl;
    }*/
  const_iterator i1=begin1();
  const_iterator i2=begin2();
  for(;i1!=end1();++i1,++i2){
    fresult<<std::fixed<<std::setprecision(prec)<<(*i1).X()<<"\t"<<(*i1).Y()<<"\t"<<(*i2).X()<<"\t"<<(*i2).Y()<<std::endl;
  }
  return true;
}


template<class T> const bool ResultMEC<T>::EcrireTexte(std::ofstream & fresult) const{
  if(_pts1.size()!=_pts2.size()){return false;}
  if(_pts1.empty()){return true;}
  /*  for(unsigned int n=0;n<_pts1.size();n++){
    fresult<<_pts1[n].X()<<"\t"<<_pts1[n].Y()<<"\t"<<_pts2[n].X()<<"\t"<<_pts2[n].Y()<<std::endl;
    }*/
  const_iterator i1=begin1();
  const_iterator i2=begin2();
  for(;i1!=end1();++i1,++i2){
    //fresult<<std::fixed<<std::setprecision(4)<<(*i1).X()<<"\t"<<(*i1).Y()<<"\t"<<(*i2).X()<<"\t"<<(*i2).Y()<<std::endl;
    fresult<<(*i1).X()<<"\t"<<(*i1).Y()<<"\t"<<(*i2).X()<<"\t"<<(*i2).Y()<<std::endl;
  }
  return true;
}


template<class T> template<class U> void ResultMEC<T>::Ecrire(std::ofstream & fresult, const T2D<U> & pt1, const T2D<U> & pt2){
#ifndef ES_RESULT_BINAIRE
  EcrireTexte(fresult,pt1,pt2);
#else
  EcrireBinaire(fresult,pt1,pt2);
#endif
}

template<class T> template<class U> void ResultMEC<T>::EcrireBinaire(std::ofstream & fresult, const T2D<U> & pt1, const T2D<U> & pt2){
  fresult.write((const char*)&(pt1.X()),sizeof(float));
  fresult.write((const char*)&(pt1.Y()),sizeof(float));
  fresult.write((const char*)&(pt2.X()),sizeof(float));
  fresult.write((const char*)&(pt2.Y()),sizeof(float));
}

template<class T> template<class U> void ResultMEC<T>::EcrireTexte(std::ofstream & fresult, const T2D<U> & pt1, const T2D<U> & pt2){
  fresult<<pt1.X()<<"\t"<<pt1.Y()<<"\t"<<pt2.X()<<"\t"<<pt2.Y()<<std::endl;
}
template<class T> template<class U> void ResultMEC<T>::EcrireTexte(std::ofstream & fresult, const T2D<U> & pt1, const T2D<U> & pt2, const int prec){
  if(prec<=0){ResultMEC<T>::EcrireTexte(fresult,pt1,pt2); return;}
  fresult<<std::fixed<<std::setprecision(prec)<<pt1.X()<<"\t"<<pt1.Y()<<"\t"<<pt2.X()<<"\t"<<pt2.Y()<<std::endl;
}

template<class T> template<class U> void ResultMEC<T>::Lire(std::ifstream & fresult, T2D<U> & pt1, T2D<U> & pt2){
#ifndef ES_RESULT_BINAIRE
  LireTexte(fresult,pt1,pt2);
#else
  LireBinaire(fresult,pt1,pt2);
#endif
}

template<class T> template<class U> void ResultMEC<T>::LireBinaire(std::ifstream & fresult, T2D<U> & pt1, T2D<U> & pt2){
  fresult.read((char*)&(pt1.x),sizeof(float));
  fresult.read((char*)&(pt1.y),sizeof(float));
  fresult.read((char*)&(pt2.x),sizeof(float));
  fresult.read((char*)&(pt2.y),sizeof(float));
}
template<class T> template<class U> void ResultMEC<T>::LireTexte(std::ifstream & fresult, T2D<U> & pt1, T2D<U> & pt2){
  fresult>>pt1.x()>>pt1.y()>>pt2.x()>>pt2.y();
}







template<class T> const bool ResultMEC<T>::LirePk0(const stringpp & adresse1, const stringpp & adresse2){
  std::ifstream f1(adresse1.c_str(),std::ios::in);
  if(!f1.good()){std::cout<<"Fichier "<<adresse1<<" introuvable\n"; return false;}
  std::ifstream f2(adresse2.c_str(),std::ios::in);
  if(!f2.good()){std::cout<<"Fichier "<<adresse2<<" introuvable\n"; return false;}
  return LirePk0(f1,f2);
}

template<class T> const bool ResultMEC<T>::LirePk0(std::ifstream & f1, std::ifstream & f2){
  if(!f1.good()){return false;}
  if(!f2.good()){return false;}
  int i; T pttmp;
  std::map<int,T> mpts1,mpts2;
  while(true){
    f1>>i;
    if(f1.eof())break;
    f1>>pttmp.x; f1>>pttmp.y;
    if(f1.eof())break;
    mpts1[i]=pttmp;
  }
  while(true){
    f2>>i;
    if(f2.eof())break;
    f2>>pttmp.x; f2>>pttmp.y;
    if(f2.eof())break;
    mpts2[i]=pttmp;
  }
  typename std::map<int,T>::const_iterator it1,it2;
  for(it1=mpts1.begin();it1!=mpts1.end();++it1){
    it2=mpts2.find((*it1).first);
    if(it2==mpts2.end()){continue;}
    _pts1.push_back((*it1).second);
    _pts2.push_back((*it2).second);
  }
}

template<class T> const bool ResultMEC<T>::EcrirePk0(const stringpp & adresse1, const stringpp & adresse2, const int imin) const{
  std::ofstream f1(adresse1.c_str(),std::ios::out);
  std::ofstream f2(adresse2.c_str(),std::ios::out);
  return EcrirePk0(f1,f2,imin);
}
template<class T> const bool ResultMEC<T>::EcrirePK0(std::ofstream & f1, std::ofstream & f2, const int imin) const{
  if(_pts1.size()!=_pts2.size()){return false;}
  if(_pts1.empty()){return true;}
  iterator i1=begin1();
  iterator i2=begin2();
  int id=imin;
  for(;i1!=end1();++i1,++i2,++id){
    f1<<id<<"\t"<<(*i1).X()<<"\t"<<(*i1).Y()<<std::endl;
    f2<<id<<"\t"<<(*i2).X()<<"\t"<<(*i2).Y()<<std::endl;
  }
}







//Renvoie l'emprise des points dans la premiere image
template<class T> const TRectangle2D<T> ResultMEC<T>::Emprise1() const{
  TRectangle2D<T> emp((*(_pts1.begin())));
  for(const_iterator i=begin1();i!=end1();++i){
    if((*i).x<emp.p1.x){emp.p1.x=(*i).x;}
    if((*i).y<emp.p1.y){emp.p1.y=(*i).y;}
    if((*i).x>emp.p2.x){emp.p2.x=(*i).x;}
    if((*i).y>emp.p2.y){emp.p2.y=(*i).y;}
  }
  return emp;
}

//Renvoie l'emprise des points dans la seconde image
template<class T> const TRectangle2D<T> ResultMEC<T>::Emprise2() const{
  TRectangle2D<T> emp((*(_pts1.begin())));
  for(const_iterator i=begin2();i!=end2();++i){
    if((*i).x<emp.p1.x){emp.p1.x=(*i).x;}
    if((*i).y<emp.p1.y){emp.p1.y=(*i).y;}
    if((*i).x>emp.p2.x){emp.p2.x=(*i).x;}
    if((*i).y>emp.p2.y){emp.p2.y=(*i).y;}
  }
  return emp;
}

#ifdef _AVECMATIS_
//Ne conserve que les points du (des) masque(s)
template<class T> void ResultMEC<T>::FiltreMasque1(const TTImage<bool> & masque){
  std::vector<T> tmp1,tmp2;
  /*  for(unsigned int n=0;n<_pts1.size();n++){
    if(masque((int)_pts1[n].X(),(int)_pts1[n].Y())!=0){tmp1.push_back(_pts1[n]);tmp2.push_back(_pts2[n]);}
    }*/
  const_iterator i1=begin1();
  const_iterator i2=begin2();
  for(;i1!=end1();++i1,++i2){
    if(masque((int)*i1.X(),(int)*i1.Y())!=0){tmp1.push_back(*i1);tmp2.push_back(*i2);}
  }
  _pts1=tmp1;_pts2=tmp2;  
}

template<class T> void ResultMEC<T>::FiltreMasque2(const TTImage<bool> & masque){
  std::vector<T> tmp1,tmp2;
  /*  for(unsigned int n=0;n<_pts1.size();n++){
    if(masque((int)_pts2[n].X(),(int)_pts2[n].Y())!=0){tmp1.push_back(_pts1[n]);tmp2.push_back(_pts2[n]);}
    }*/
  const_iterator i1=begin1();
  const_iterator i2=begin2();
  for(;i1!=end1();++i1,++i2){
    if(masque((int)*i2.X(),(int)*i2.Y())!=0){tmp1.push_back(*i1);tmp2.push_back(*i2);}
  }
  _pts1=tmp1;_pts2=tmp2;  
}

template<class T> void ResultMEC<T>::FiltreMasque(const TTImage<bool> & masque1, const TTImage<bool> & masque2){
  std::vector<T> tmp1,tmp2;
  /*  for(unsigned int n=0;n<_pts1.size();n++){
    if(masque1((int)_pts1[n].X(),(int)_pts1[n].Y())!=0 && masque2((int)_pts2[n].X(),(int)_pts2[n].Y())!=0){tmp1.push_back(_pts1[n]);tmp2.push_back(_pts2[n]);}
    }*/
  const_iterator i1=begin1();
  const_iterator i2=begin2();
  for(;i1!=end1();++i1,++i2){
    if(masque1((int)(*i1).X(),(int)(*i1).Y())!=0 && masque2((int)(*i2).X(),(int)(*i2).Y())!=0){tmp1.push_back(*i1);tmp2.push_back(*i2);}
  }
  _pts1=tmp1;_pts2=tmp2;  
}

#endif















#ifdef _AVECMATIS_

//Classe heritee capable d'egalement charger les parametres de prise de vue (fichier .ori, ajout de fonction liee a la distortion et au changement de repere (pk1...)

template<class T>
class ResultMECori:public ResultMEC<T>{
 public:
  ResultMECori():ResultMEC<T>(){}
  ResultMECori(const stringpp & adresseresult):ResultMEC<T>(adresseresult){}
  ResultMECori(const stringpp & adresseresult, const stringpp & nom1, const stringpp & nom2):ResultMEC<T>(adresseresult,nom1,nom2){}
  ResultMECori(const std::vector<T> & pts1, const std::vector<T> & pts2):ResultMEC<T>(pts1,pts2){}
  ResultMECori(const stringpp & nom1, const stringpp & nom2):ResultMEC<T>(nom1,nom2){}
  ResultMECori(const std::vector<T> & pts1, const std::vector<T> & pts2, const stringpp & nom1, const stringpp & nom2):ResultMEC<T>(nom1,pts1,nom2,pts2){}
  ResultMECori(const stringpp & adresseresult, const Orientation & ori1, const Orientation & ori2):ResultMEC<T>(adresseresult),_ori1(ori1),_ori2(ori2){}
  ResultMECori(const stringpp & adresseresult, const stringpp & nom1, const stringpp & nom2, const Orientation & ori1, const Orientation & ori2):ResultMEC<T>(adresseresult,nom1,nom2),_ori1(ori1),_ori2(ori2){}
  ResultMECori(const std::vector<T> & pts1, const std::vector<T> & pts2, const Orientation & ori1, const Orientation & ori2):ResultMEC<T>(pts1,pts2),_ori1(ori1),_ori2(ori2){}
  ResultMECori(const stringpp & nom1, const stringpp & nom2, const Orientation & ori1, const Orientation & ori2):ResultMEC<T>(nom1,nom2),_ori1(ori1),_ori2(ori2){}
  ResultMECori(const std::vector<T> & pts1, const std::vector<T> & pts2, const stringpp & nom1, const stringpp & nom2, const Orientation & ori1, const Orientation & ori2):ResultMEC<T>(nom1,pts1,nom2,pts2),_ori1(ori1),_ori2(ori2){}
  ~ResultMECori(){};
  
  const Orientation ori1() const{return _ori1;}
  Orientation & ori1(){return _ori1;}
  const Orientation ori2() const{return _ori2;}
  Orientation & ori2(){return _ori2;}
  
  void Distortion_Ground2Photo();
  void Distortion_Photo2Ground();
  void RepereImage2ReperePhotogrammetrique();
  void ReperePhotogrammetrique2RepereImage();
  
 protected:
  Orientation _ori1,_ori2;
};

#endif

template<class T>
class ResultMECindice:public ResultMEC<T>{
 public:
  ResultMECindice():ResultMEC<T>(){}
  ResultMECindice(const stringpp & adresseresult):ResultMEC<T>(adresseresult){initindice();}
  ResultMECindice(const stringpp & adresseresult, const stringpp & nom1, const stringpp & nom2):ResultMEC<T>(adresseresult,nom1,nom2){initindice();}
  ResultMECindice(const std::vector<T> & pts1, const std::vector<T> & pts2):ResultMEC<T>(pts1,pts2){initindice();}
  ResultMECindice(const stringpp & nom1, const stringpp & nom2):ResultMEC<T>(nom1,nom2){}
  ResultMECindice(const std::vector<T> & pts1, const std::vector<T> & pts2, const stringpp & nom1, const stringpp & nom2):ResultMEC<T>(nom1,pts1,nom2,pts2){initindice();}
  ResultMECindice(const std::vector<int> & indices, const std::vector<T> & pts1, const std::vector<T> & pts2):ResultMEC<T>(pts1,pts2),_indices(indices){}
  ResultMECindice(const std::vector<int> & indices, const std::vector<T> & pts1, const std::vector<T> & pts2, const stringpp & nom1, const stringpp & nom2):ResultMEC<T>(nom1,pts1,nom2,pts2),_indices(indices){}
 protected:
  std::vector<int> _indices;
  void initindice(){for(unsigned int n=1;n<ResultMEC<T>::_pts1.size()+1;n++){_indices.push_back(n);}}
};


#endif
