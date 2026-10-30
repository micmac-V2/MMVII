#ifndef __TRECTANGLE2D_HPP__
#define __TRECTANGLE2D_HPP__

#include <fstream>
#include "t2d.hpp"

template<class T>
class TRectangle2D{
public :
  //Constructeurs
  TRectangle2D():_p1(),_p2() {} ;
  TRectangle2D(const T2D<T> & p):_p1(p), _p2(p) {/*_p2.x()++;_p2.y()++;*/}
  TRectangle2D(const T x, const T y):_p1(x,y), _p2(x,y) {/*_p2.x()++;_p2.y()++;*/}
  TRectangle2D(const T x1, const T y1, const T x2, const T y2):_p1(std::min(x1,x2),std::min(y1,y2)),_p2(std::max(x1,x2),std::max(y1,y2)) {} 
  TRectangle2D(const TRectangle2D & rec):_p1(rec._p1),_p2(rec._p2) {} 
  TRectangle2D(const T2D<T> & a, const T2D<T> & b):_p1(std::min(a.x(),b.x()),std::min(a.y(),b.y())),_p2(std::max(a.x(),b.x()),std::max(a.y(),b.y())) {}
  template<class U> //U=T2D<T> ou T3D>T>
  TRectangle2D(const std::vector<U> & v){
    typename std::vector<U>::const_iterator it=v.begin();
    _p1=T2D<T>((*it).x(),(*it).y());
    _p2=_p1;
    for(;it!=v.end();++it){Ajoute((*it).x(),(*it).y());}
  }
  
  //Rectangle contient point ?
  template<class U> const bool contient (const T2D<U> & p) const{
    return (p.X()>=_p1.X() && p.Y()>=_p1.Y() && p.X()<=_p2.X() && p.Y()<=_p2.Y()) ;
  }
  const bool contient (const T2D<T> & p) const ; 
  template<class U> const bool contient (const U x, const U y) const{
    return ((double)x>=_p1.X() && (double)y>=_p1.Y() && (double)x<=_p2.X() && (double)y<=_p2.Y()) ;
  }
  
  
  //Union de rectangles
  const TRectangle2D<T> operator | (const TRectangle2D<T> & r) const {
    return TRectangle2D<T>(std::min(_p1.x(),r._p1.x()),std::min(_p1.y(),r._p1.y()),std::max(_p2.x(),r._p2.x()),std::max(_p2.y(),r._p2.y())) ;
  }
  TRectangle2D<T> & operator |= (const TRectangle2D<T> & r) {
    if (r._p1.x()<_p1.x()) _p1.x()=r._p1.x() ;
    if (r._p1.y()<_p1.y()) _p1.y()=r._p1.y() ;
    if (r._p2.x()>_p2.x()) _p2.x()=r._p2.x() ;
    if (r._p2.y()>_p2.y()) _p2.y()=r._p2.y() ;
    return (*this) ;
  }
  // Augmentation d'un rectangle pour contenir un point
  const TRectangle2D<T> operator | (const T2D<T> & p) const {
    return TRectangle2D<T>( std::min(_p1.x(),p.x()),std::min(_p1.y(),p.y()),std::max(_p2.x(),p.x()+1),std::max(_p2.y(),p.y()+1)) ;
  }
  void Ajoute(const T2D<T> & p)  {
    if(p.x()<_p1.x()){_p1.x()=p.x();}
    if(p.y()<_p1.y()){_p1.y()=p.y();}
    if(p.x()>_p2.x()){_p2.x()=p.x();}
    if(p.y()>_p2.y()){_p2.y()=p.y();}
  }
  void Ajoute(const T x, const T y) {
    if(x<_p1.x()){_p1.x()=x;}
    if(y<_p1.y()){_p1.y()=y;}
    if(x>_p2.x()){_p2.x()=x;}
    if(y>_p2.y()){_p2.y()=y;}
  }

  // Changement en passant en parametre un TRectangle2D<T> au lieu de <int>
  const bool inclusDans (const TRectangle2D < T > & src) {
    return _p1.x()>=src._p1.x() && _p1.y()>=src._p1.y() && _p2.x()<=src._p2.x() && _p2.y()<=src._p2.y() ;
  }
  
  // Intersection de rectangles
  const TRectangle2D<T> operator & (const TRectangle2D<T> & r) {
    T2D<T> tmp_p1 ( std::max(_p1.x(), r._p1.x()), std::max(_p1.y(), r._p1.y())) ;
    T2D<T> tmp_p2 ( std::min(_p2.x(), r._p2.x()), std::min(_p2.y(), r._p2.y())) ;
    if (tmp_p2.x()<tmp_p1.x() || tmp_p2.y()<tmp_p1.y()) 
      return TRectangle2D<T>(0,0,0,0);
    return TRectangle2D<T>(tmp_p1,tmp_p2);
  }
  
  // vide ?
  const bool vide () const {return (_p1.x()>=_p2.x() || _p1.y() >= _p2.y());}
  const bool empty() const{return vide();}
  
  //-- Dilatation
  const TRectangle2D<T> dilater (const T v) const {
    return TRectangle2D<T> (_p1.x()-v,_p1.y()-v,_p2.x()+v,_p2.y()+v) ;
  }

  void dilaterPlace (const T v) {
    _p1.x() -= v ;
    _p1.y() -= v ;
    _p2.x() += v ;
    _p2.y() += v ;
  }
  
  
  const TRectangle2D<T> Translater (const T2D<T> & t) const{
    return TRectangle2D<T>(_p1+t,_p2+t);
  }
  void translaterPlace(const T2D<T> & t){
    _p1+=t;
    _p2+=t;
  }
  
  //-----------------------------------
  // IO
  //-----------------------------------
  void Charger(std::ifstream & i){
    T2D<T> _p1,_p2;
    i>>_p1.x();
    i>>_p1.y();
    i>>_p2.x();
    i>>_p2.y();
    *this=TRectangle2D<T>(_p1,_p2);
  }
  
  void Charger(char * nom)
  {
    std::ifstream i(nom,std::ios::in) ;
    Charger(i);
  }

  const T2D<T> & p1() const{return _p1;}
  T2D<T> & p1() {return _p1;}
  const T2D<T> & p2() const{return _p2;}
  T2D<T> & p2() {return _p2;}
  //Ajout perso pour parcours sommets...
  const T2D<T> operator[] (const int n) const{
    if(n==0){
      return _p1;
    }else if(n==1){
      return T2D<int>(_p1.x(),_p2.y());
    }else if(n==2){
      return _p2;
    }else if(n==3){
      return T2D<int>(_p2.x(),_p1.y());
    }
    if(n>3){return T2D<int>(_p2.x(),_p1.y());}
    if(n<0){return _p1;}
    return _p1;
  }
  
  void Afficher(std::ostream & fout) const{
    fout<<_p1<<" "<<_p2<<std::endl;
  }


protected:

  T2D<T> _p1, _p2 ;
  
} ;

template <>
inline const bool TRectangle2D < int >::contient (const T2D < int > & p) const {
  return (p.x()>=_p1.x() && p.y()>=_p1.y() && p.x()<_p2.x() && p.y()<_p2.y()) ;
}

template <>
inline const bool TRectangle2D < double >::contient (const T2D < double >  & p) const {
  return (p.x()>=_p1.x() && p.y()>=_p1.y() && p.x()<=_p2.x() && p.y()<=_p2.y()) ;
}



template < class T, class U >
void caster (TRectangle2D < T > &t, TRectangle2D < U > & u) {
  t._p1= (T2D<T>) u._p1 ;
  t._p2= (T2D<T>) u._p2 ;
}


template <class T>
std::ostream& operator << ( std::ostream &s, const TRectangle2D<T> &r ) {
  s<<"("<<r._p1<<","<<r._p2<<")" ;
  return (s) ;
}

//=======================================================
//Parcourir un rectangle.
//Pour a et b entier. Les coordonnees (a,b) parcourent toutes les positions possibles dans le rectangle r.
//Le parcours est fait par lignes.
//-------------------------------------------------------
#ifndef __TRECTANGLE_HPP__
#define BOUCLERECTANGLE(a,b,r) for((b)=r._p1.y();(b)<r._p2.y();++(b)) for((a)=r._p1.x();(a)<r._p2.x();++(a)) 
#endif

#endif
