#ifndef _T2D_
#define _T2D_

#include <iostream>
#include <cmath>
#include <vector>

#ifdef _T3D_
template<class T> class T3D;
#endif

template<class T> 
class T2D{
protected:
  T _x;
  T _y;
  
public:
  //Constructeurs et destructeurs
  T2D():_x(0),_y(0){};
  T2D(const T v):_x(v),_y(v){};
  T2D(const T X, const T Y):_x(X),_y(Y){};
  T2D(const T2D & src):_x(src.x()),_y(src.y()){}
  template<class V> T2D(const V X,const V Y):_x((T)X),_y((T)Y){};
  template<class V> T2D(const T2D<V> & src):_x((T)src.x()),_y((T)src.y()){}
  template<class V> T2D(const std::vector<V> & src):_x((T)src[0]),_y((T)src[1]){}
  
  virtual ~T2D(){};

  const unsigned int dim() const{return 2;}
  const unsigned int size() const{return 2;}

  //typedef
  typedef T type;

  //Accesseurs
  inline const T x() const{return _x;} 
  inline T & x() {return _x;}
  inline const T y() const{return _y;}
  inline T & y() {return _y;}
  inline const T p1() const{return _x;} 
  inline T & p1() {return _x;}
  inline const T p2() const{return _y;}
  inline T & p2() {return _y;}
  inline const T first() const{return _x;} 
  inline T & first() {return _x;}
  inline const T second() const{return _y;}
  inline T & second() {return _y;}
  inline const T real() const{return _x;} 
  inline T & real() {return _x;}
  inline const T imag() const{return _y;}
  inline T & imag() {return _y;}

  const T operator [] (const int n) const {
    if (n==0) return x() ;
    if (n==1) return y() ;
    return x() ;
  }
  T & operator [] (const int n) {
    if (n==0) return x() ;
    if (n==1) return y() ;
    return x() ;
  }

  //Conservé pour compatibilité
  template<class V> inline void setx(const V v) { x()=(T)v;}
  template<class V> inline void sety(const V v) { y()=(T)v;}
  inline void set (const T a, const T b) {x()=a; y()=b;}
  template<class V> inline void set (const V a, const V b) {x()=a; y()=b;}
  template<class V> inline void set (const std::vector<V> & v) {x()=v[0]; y()=v[1];}

  
  //Pour travailler en double
  inline const double X() const{return (double)x();}
  inline const double Y() const{return (double)y();}
  inline const T2D<double> todouble() const{return T2D<double>(X(),Y());}
  inline const std::vector<T> tovector() const{std::vector<T> v(2); v[0]=x(); v[1]=y(); return v;}

  //Pour travailler avec autre chose
  template<class V> void cast(T2D<V> & p) const{p.set(x(),y());}
  template<class V> void cast(std::vector<V> & p) const{p.resize(2); p[0]=x(); p[1]=y();}
#ifdef _T3D_
  const T3D<T> castT3D() const{return T3D<T>(x(),y(),0);}
#endif

  //Operateurs
  const bool operator == (const T2D<T> & p) const { return (x()==p.x() && y()==p.y()) ; }
  const bool operator != (const T2D<T> & p) const { return (x()!=p.x() || y()!=p.y()) ; }
  const bool operator < (const T2D<T> & p) const { return ( (x()<p.x()) || ((x()==p.x())&&(y()<p.y())) ) ; }
  const bool operator > (const T2D<T> & p) const { return ( (x()>p.x()) || ((x()==p.x())&&(y()>p.y())) ) ; }

  //Autres
  const T2D<T> swapxy() const {return (T2D<T> (y(),x()));}


  //Operateurs
  const T2D operator - (const T2D & b) const {return T2D<T>(x()-b.x(),y()-b.y());}
  const T2D operator - () const {return T2D<T>(-x(),-y());}
  const T2D operator + (const T2D &b) const {return T2D(x()+b.x(),y()+b.y());}
  void operator += (const T2D &b) {x()+=b.x(); y()+=b.y();} 
  const T2D operator / (const T t) const {return T2D(x()/t,y()/t);}
  const T2D operator * (const T t) const {return T2D(x()*t,y()*t);}
  T2D & operator /= (const T t) {x()/=t; y()/=t; return(*this);}  
  T2D & operator *= (const T t) {x()*=t; y()*=t; return (*this);}

  //Produit scalaire
  const double operator * (const T2D & p) const {return X()*p.X()+Y()*p.Y();}
  // Produit vectoriel (déterminant)
  const double operator ^ (const T2D & p) const {return X()*p.Y()-Y()*p.X();}

  /// Norme
  const double norme() const {return sqrt(X()*X()+Y()*Y());}
  const double norme2() const {return (X()*X()+Y()*Y());}

  T2D<T> & normaliser () {
    double d=norme() ;
    if (d){
      x()=(T)((double)x()/d);
      y()=(T)((double)y()/d);
    }
    return (*this);
  }
  void normalise(){normaliser();}

  template<class V> const double dist2(const T2D<V> & p) const{
    double d1=p.X()-X();
    double d2=p.Y()-Y();
    return d1*d1+d2*d2;
  }
  template<class V> const double dist(const T2D<V> & p) const{
    return sqrt(dist2(p));
  }

  template <class U, class V> 
  static inline const double dist(const T2D<U> & p1, const T2D<V> & p2){
    return p1.dist(p2);
  }
  template <class U, class V> 
  static inline const double dist2(const T2D<U> & p1, const T2D<V> & p2){
    return p1.dist2(p2);
  }


};

template <class T, class U> 
inline T2D<T> & operator -= (T2D<T> &a, const T2D<U> &b) {
  a.x() -= b.x();
  a.y() -= b.y();
  return a;
}

template <class T , class V>
inline T2D<T> operator * (const V t,const T2D<T> & p) {return p*t;}

template <class T, class V>
inline const double dist(const T2D<T> & p1, const T2D<V> & p2){
  return p1.dist(p2);
}

template <class T>
inline std::ostream& operator << ( std::ostream & os, const T2D<T> & p ) {
  return (os <<"("<<p.x()<<","<<p.y()<<")");
}




  
 
#endif
