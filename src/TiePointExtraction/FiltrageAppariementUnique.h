#ifndef _FILTRAGEAPPARIEMENTUNIQUE_HPP_
#define _FILTRAGEAPPARIEMENTUNIQUE_HPP_


#include <set>
#include <vector>


//=============================================================
//Classe FiltreAppariementUnique
// Filtrage des appariements en verifiant par aller retour que deux points distincts ne sont pas apparies a un meme point
//=============================================================
class FiltreAppariementUnique{
 protected:
  bool _doubloninterdit;
 public:
  //Constructeur / Destructeur
 FiltreAppariementUnique(const bool doubloninterdit=false):_doubloninterdit(doubloninterdit){};
  virtual ~FiltreAppariementUnique(){}

  //Accesseur
  const bool doubloninterdit() const{return _doubloninterdit;}
  bool & doubloninterdit(){return _doubloninterdit;}
  
  //Fonctions
  template<class T>
    void filtre(std::vector<T> & pts1, std::vector<T> & pts2) const{
    FiltreAppariementUnique::FiltrageAppariementUnique(pts1,pts2,_doubloninterdit);
  }
  template<class T>
    void operator() (std::vector<T> & pts1, std::vector<T> & pts2) const{
    filtre(pts1,pts2);
  }


  //Fonctions statics

  //------------------------------------------------------------------
  //Effectue un filtrage des points faux en verifiant l'unicite de l'appariement
  //------------------------------------------------------------------
  template<class T>
  static void FiltrageAppariementUnique(std::vector<T> & pts1, std::vector<T> & pts2, const bool doubloninterdit=false){
    if(pts1.empty() || pts2.empty()){return;}
    std::vector<T> ptstmp1,ptstmp2;
    T2D<double> xy2,xy1;
    std::set<int> avirer;
    for(unsigned int n=0;n<pts1.size();n++){
      xy2=T2D<double>(pts2[n].X(),pts2[n].Y());
      typename std::vector<T>::const_iterator it2=pts2.begin()+n+1;
      xy1=T2D<double>(pts1[n].X(),pts1[n].Y());
      typename std::vector<T>::const_iterator it1=pts1.begin()+n+1;
      for(unsigned int m=n+1;m<pts2.size();m++,++it1,++it2){
	if((*it2).X()==xy2.x() && (*it2).Y()==xy2.y()){
	  if(doubloninterdit){
	    avirer.insert(m); 
	    continue;
	  }
	  if(pts1[m].X()!=pts1[n].X() || pts1[m].Y()!=pts1[n].Y()){
	    avirer.insert(n);
	    avirer.insert(m);
	  }
 	}
	if((*it1).X()==xy1.x() && (*it1).Y()==xy1.y()){
	  if(doubloninterdit){
	    avirer.insert(m); 
	    continue;
	  }
	  if(pts2[m].X()!=pts2[n].X() || pts2[m].Y()!=pts2[n].Y()){
	    avirer.insert(n);
	    avirer.insert(m);
	  }
	}
      }
    }
    for(unsigned int n=0;n<pts1.size();n++){
      if(avirer.find(n)!=avirer.end()){continue;}
      ptstmp1.push_back(pts1[n]);
      ptstmp2.push_back(pts2[n]);
    }
    pts1=ptstmp1;
    pts2=ptstmp2;
  }


  //------------------------------------------------------------------
  //Supprime simplement les doublons
  //------------------------------------------------------------------
  template<class T>
  static void FiltrageDoublons(std::vector<T> & pts1, std::vector<T> & pts2){
    if(pts1.empty() || pts2.empty()){return;}
    std::vector<T> ptstmp1,ptstmp2;
    T2D<double> xy2,xy1;
    std::set<int> avirer;
    for(unsigned int n=0;n<pts1.size();n++){
      xy2=T2D<double>(pts2[n].X(),pts2[n].Y());
      xy1=T2D<double>(pts1[n].X(),pts1[n].Y());
      typename std::vector<T>::const_iterator it2=pts2.begin()+n+1;
      typename std::vector<T>::const_iterator it1=pts1.begin()+n+1;
      for(unsigned int m=n+1;m<pts2.size();m++,++it1,++it2){
	if((*it2).X()==xy2.x() && (*it2).Y()==xy2.y()){avirer.insert(m);}
	if((*it1).X()==xy1.x() && (*it1).Y()==xy1.y()){avirer.insert(m);}
      }
    }
    for(unsigned int n=0;n<pts1.size();n++){
      if(avirer.find(n)!=avirer.end()){continue;}
      ptstmp1.push_back(pts1[n]);
      ptstmp2.push_back(pts2[n]);
    }
    pts1=ptstmp1;
    pts2=ptstmp2;
  }


};


#endif
