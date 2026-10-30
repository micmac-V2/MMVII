#ifndef _EMPRISESCLICHECLICHE_HPP_
#define _EMPRISESCLICHECLICHE_HPP_

#include "TPolygone2D.hpp"


class EmprisesClicheCliche{
public:
  EmprisesClicheCliche(){}
  EmprisesClicheCliche(const std::string & adresse_fichier_emprises){Lire(adresse_fichier_emprises);}
  virtual ~EmprisesClicheCliche(){}

  //Accesseurs
  const TPolygone2D<double> & emprise_2_dans_1() const{return _emprise_2_dans_1;}
  TPolygone2D<double> & emprise_2_dans_1(){return _emprise_2_dans_1;}
  const TPolygone2D<double> & emprise_1_dans_2() const{return _emprise_1_dans_2;}
  TPolygone2D<double> & emprise_1_dans_2(){return _emprise_1_dans_2;}

  //Lecture
  const bool Lire(const std::string & adresse_fichier_emprises){
    return LectureEmprisesTexte(adresse_fichier_emprises,_emprise_2_dans_1,_emprise_1_dans_2);
  }

  //Methodes statiques
  template<class T> static const bool LectureEmprisesTexte(const std::string & adresse_fichier_emprises, TPolygone2D<T> & emprise_2_dans_1, TPolygone2D<T> & emprise_1_dans_2);
  
protected:
  TPolygone2D<double> _emprise_2_dans_1;
  TPolygone2D<double> _emprise_1_dans_2;
};



template<class T>
const bool EmprisesClicheCliche::LectureEmprisesTexte(const std::string & adresse_fichier_emprises, TPolygone2D<T> & emprise_2_dans_1, TPolygone2D<T> & emprise_1_dans_2){
  std::ifstream femp(adresse_fichier_emprises.c_str(),std::ios::in);
  if(!femp.good()){
    std::cout<<"Impossible de lire le fichier '"<<adresse_fichier_emprises<<"'"<<std::endl;
    return false;
  }
  unsigned int nbpts;
  femp>>nbpts;
  std::vector<T2D<T> > emptmp(nbpts);
  for(typename std::vector<T2D<T> >::iterator i=emptmp.begin();i!=emptmp.end();++i){
    femp>>(*i).x()>>(*i).y();
  }
  emprise_2_dans_1=TPolygone2D<T>(emptmp);
  emptmp.clear();
  femp>>nbpts;
  emptmp.resize(nbpts);
  for(typename std::vector<T2D<T> >::iterator i=emptmp.begin();i!=emptmp.end();++i){
    femp>>(*i).x()>>(*i).y();
  }
  emprise_1_dans_2=TPolygone2D<T>(emptmp);
  
  return true;
}
 

#endif
