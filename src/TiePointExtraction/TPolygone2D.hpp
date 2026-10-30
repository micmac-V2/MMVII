#ifndef TPOLYGONE2D_HPP
#define TPOLYGONE2D_HPP

#include <iostream>
#include <vector>
#include "trectangle.hpp"

template<class T>
class TPolygone2D {
protected:
  std::vector<T2D<T> > _points;
  
public:
  TPolygone2D() {;}
  TPolygone2D(const std::vector<T2D<T> > & pts):_points(pts) {;}
  virtual ~TPolygone2D(){}
  
  // Ajout d'un point
  void AjoutePoint(const T2D<T> & M) { _points.push_back(M);}
  void AjoutePoint(const double x, const double y) { _points.push_back(T2D<T>(x,y));}
  
  const unsigned int nb_points() const{ return _points.size();}
  const unsigned int size() const{ return _points.size();}
  const T2D<T> Pt(const unsigned int num) const{ if(num<_points.size()) return _points[num]; else return T2D<T>(0,0);}

  //Accesseurs
  const std::vector<T2D<T> > & points() const{return _points;}
  std::vector<T2D<T> > & points() {return _points;}

  
  // Suppression de tous les points
  void clear() { _points.clear();}
  
  // Suppression du dernier point
  void SupprimeDernierPoint() { _points.pop_back();}
  
  // Fermeture du polygone
  void fermer(){
    if(est_ferme()){return;}
    //if(size()<3){return;}
    _points.push_back(*(_points.begin()));
  }
  
  // Test si un polygone est ferme
  bool est_ferme() const{ 
    if(_points.size()<3) return false; 
    if(_points[0]==_points[size()-1]) return true; 
    return false;
  }
  
  // Retourne le rectangle englobant du polygone
  const TRectangle2D<T> Emprise() const;

  // Retourne le polygone convexe englobant le polygone
  const TPolygone2D convexe() const;
  
  // Indique si un point est contenu dans le polygone
  const bool contient(const T2D<T> & P) const;
  const bool contient(const T x, const T y) const {return contient(T2D<T>(x,y));}
  
  // Regarde si deux polygones s'intersectent
  const bool intersecte(const TPolygone2D & polygone) const;
  
  // Renvoie le barycentre du polygone
  const T2D<double> barycentre() const;
  
  // Expansion du polygone
  const TPolygone2D dilate(const double d) const;
  
  // Surface du polygone
  const double aire() const;
  const double surface() const{return aire();}

  // Longueur du polygone
  const double perimetre() const;
  const double longueur() const{return perimetre();}

  // Retourne true si le polygone est decrit dans le sens trigonometrique, false sinon
  const bool sens_trigonometrique() const;
  //change le sens...
  void inverse_sens();
  const TPolygone2D get_sens_inverse() const;
  void set_sens(const bool trigo);
  
};


//-----------------------------------------------------------------------------
// Retourne le rectangle englobant du polygone
//-----------------------------------------------------------------------------
template<class T>
const TRectangle2D<T> TPolygone2D<T>::Emprise() const{
  if(size()<1) return TRectangle2D<T>(0,0,0,0);
  
  double xmin, xmax, ymin, ymax;
  T2D<T> ptmp(_points[0]);
  xmin=xmax=ptmp.x();
  ymin=ymax=ptmp.y();
  
  for(unsigned int i=1;i<size();i++){
    ptmp=_points[i];
    xmin=std::min(xmin,ptmp.x());
    ymin=std::min(ymin,ptmp.y());
    xmax=std::max(xmax,ptmp.x());
    ymax=std::max(ymax,ptmp.y());
  }
  return TRectangle2D<T>(xmin,ymin,xmax,ymax);
}

//-----------------------------------------------------------------------------
// Retourne le polygone convexe englobant le polygone
//-----------------------------------------------------------------------------
template<class T>
const TPolygone2D<T> TPolygone2D<T>::convexe() const{
  //TODO!!!!!
  throw("Fonction a ecrire");
  std::cout<<"Fonction a ecrire...";
  return TPolygone2D<T>();
}

//-----------------------------------------------------------------------------
// Indique si un point est contenu dans le polygone
//-----------------------------------------------------------------------------
template<class T>
const bool TPolygone2D<T>::contient(const T2D<T> & pt) const{
    if(size()<3) return false;
    TRectangle2D<T> emprise=Emprise();
    if(!emprise.contient(pt)) return false;

    T2D<T> M(pt.X(),pt.Y());
    double alpha,delta=0.0;
    int nb_inter=0;
    T2D<double> A,B,C;
    for(unsigned int i=0;i<_points.size()-1;i++){
      A=_points[i];
      B=_points[i+1];
      if((A.X()==M.X())&&(A.Y()==M.Y())){
	nb_inter=0;
	break;
      }
      alpha=((A.X()-A.Y())-(M.X()+delta-M.Y()))/((A.X()-A.Y())-(B.X()-B.Y()));
      if(alpha==0.0){
	delta+=0.05;
	nb_inter = 0;
	i=-1;
	continue;
      }

      if((alpha<=0.0)||(alpha>1.0))continue;
      C=A+alpha*(B-A);
      C-=M;
      if((C.X()>0.0)&&(C.Y()>0.0))
	nb_inter++;
    }
    if((nb_inter%2)==0){    // Nombre pair d'intersection
      return false;
    }else{
      return true;
    }
}

//-----------------------------------------------------------------------------
// Regarde si deux emprises s'intersectent
//-----------------------------------------------------------------------------
template<class T>
const bool TPolygone2D<T>::intersecte(const TPolygone2D<T> & pol) const{
  for(unsigned int i=0;i<pol.size();i++){
    if(contient(pol.Pt(i))){
      return true;
    }
  }
  return false;
}

//-----------------------------------------------------------------------------
// Renvoie le barycentre du polygone
//-----------------------------------------------------------------------------
template <class T>
const T2D<double> TPolygone2D<T>::barycentre() const{
  T2D<double> bary;
  for(typename std::vector<T2D<T> >::const_iterator i=_points.begin();i!=_points.end();++i){
    //bary += (_points[i]/_points.size());
    bary.x()=bary.x()+(*i).X();
    bary.y()=bary.y()+(*i).Y();
  }
  return bary/_points.size();
}


//-----------------------------------------------------------------------------
// Expansion du polygone
//-----------------------------------------------------------------------------
template<class T>
const TPolygone2D<T> TPolygone2D<T>::dilate(const double d) const{
  T2D<double> bary=barycentre();
  TPolygone2D<T> poly;
  
  for(unsigned int i=0;i<_points.size();i++) {
    T2D<double> M(_points[i]);
    poly.AjoutePoint(((dist(bary, M) + d)/dist(bary, M))*(M-bary) + bary);
  }
  return poly;
}

//-----------------------------------------------------------------------------
// Surface du polygone
//-----------------------------------------------------------------------------
template<class T>
const double TPolygone2D<T>::aire()const{
  if(size()<3){return 0.;}
  double area=0.;
  for(unsigned int i = 0;i<_points.size()-1;i++){
    area += ((_points[i+1].Y+_points[i].Y)*(_points[i+1].X-_points[i].X)*0.5);
  }
  // Attention : il faut fermer sur le premier point du polygone
  unsigned int i = _points.size()-1;
  area += ((_points[0].Y+_points[i].Y) * (_points[0].X-_points[i].X) * 0.5);
  return area;
}

//-----------------------------------------------------------------------------
// Perimetre du polygone
//-----------------------------------------------------------------------------
template<class T>
const double TPolygone2D<T>::perimetre() const{
  double d=0.;
  if(size()<2){return 0;}
  for (unsigned int i=0;i<_points.size()-1;i++){d += dist2(_points[i],_points[i+1]);}
  return sqrt(d);
}

//-----------------------------------------------------------------------------
// Retourne true si le polygone est decrit dans le sens trigonometrique, false sinon
//-----------------------------------------------------------------------------
template<class T>
const bool TPolygone2D<T>::sens_trigonometrique() const{
  if(_points.size()<3){return true; }
  unsigned int nbPt = _points.size();
  if(est_ferme()){nbPt -= 1;}
  double xmin = _points[0].X;
  unsigned int num = 0;
  for(unsigned int i = 1; i < nbPt; i++)
    if(_points[i].X < xmin) {
      xmin = _points[i].X;
      num = i;
    }

  T2D<double> A, B = _points[num], C;
  if(num == 0){
    A = _points[nbPt-1];
  }else{
    A = _points[num-1];
  } 
  if(num == nbPt-1){
    C = _points[0];
  }else{
    C = _points[num+1];
  }
  if(prodMixt(T3D<double>(A,0),T3D<double>(B,0),T3D<double>(C,0)) > 0){return true;}
  return false;
}

//-----------------------------------------------------------------------------
//change le sens...
//-----------------------------------------------------------------------------
template<class T>
void TPolygone2D<T>::inverse_sens(){
  std::vector<T2D<T> > points(_points);
  int ilast=_points.size()-1;
  for(int n=0;n<_points.size();n++){
    _points[ilast-n]=points[n];
  }
}

template<class T>
const TPolygone2D<T> TPolygone2D<T>::get_sens_inverse() const{
  std::vector<T2D<T> > points(_points.size());
  int ilast=_points.size()-1;
  for(int n=0;n<_points.size();n++){
    points[ilast-n]=_points[n];
  }
  return TPolygone2D<T>(points);
}

template<class T>
void TPolygone2D<T>::set_sens(const bool trigo){
  bool sensactuel=sens_trigonometrique();
  if(sensactuel && trigo){return;}
  if(!sensactuel && !trigo){return;}
  inverse_sens();
}


#endif
