#ifndef _ACCESIMAGE_HPP_
#define _ACCESIMAGE_HPP_
#include "MMVII_DeclareCste.h"
#include "MMVII_Tpl_Images.h"
#include "t2d.hpp"

class AccesImage{
protected:
  MMVII::cDataFileIm2D _aFileIn;
public:
  //  AccesImage(){}
  AccesImage(const std::string & adresse_image):_aFileIn(MMVII::cDataFileIm2D::Create(adresse_image,MMVII::eForceGray::No)){;}
   ~AccesImage(){}
  inline const unsigned int GetNbColonnes() const{return _aFileIn.Sz().x();}
  inline const unsigned int GetNbLignes() const{return _aFileIn.Sz().y();}
  inline const T2D<unsigned int> GetTaille() const{return T2D<unsigned int>(_aFileIn.Sz().x(),_aFileIn.Sz().y());}
  inline const unsigned int GetNbCanaux() const{return 1;}
  inline const bool EstDans(const int i, const int j) const{return (i>=0 && j>=0 && i<(int)GetNbColonnes() && j<(int)GetNbLignes());}
  const MMVII::cDataFileIm2D & operator()() const{return _aFileIn;}

};

typedef MMVII::cIm2D<MMVII::tREAL4> ImageMemoire;


#endif

/*
  cDataFileIm2D aFileIn= cDataFileIm2D::Create(mNameIn,eForceGray::No);
  std::cout<<aFileIn.Sz()<<std::endl;
  std::cout<<aFileIn.Sz().x()<<std::endl;
  std::cout<<aFileIn.Sz().y()<<std::endl;
*/
