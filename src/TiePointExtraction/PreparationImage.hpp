#ifndef _PREPARATIONIMAGE_HPP_
#define _PREPARATIONIMAGE_HPP_
#include "AccesImage.hpp"
#include "sift.hpp"
#include "Parametres_Pretraitement_Dallage.hpp"


//=================================================
//classe PreparationImage 
// Quelques fonctions de traitement afin de preparer les images :
// - crop de la dalle
// - pretraitements (choix du canal, calcul du gradient, reetalement...) parametres par la classe ParametragePretraitement
// - conversion du format externe (classe ImageMemoire) vers le format interne (VLImageBuffer)
//-------------------------------------------------


class PreparationImage{
 protected:
  ParametragePretraitement _parametres;
 public:
  PreparationImage(){}
  PreparationImage(ParametragePretraitement & parametres):_parametres(parametres){}
  virtual ~PreparationImage(){}

  const ParametragePretraitement & parametres() const {return _parametres;}
  ParametragePretraitement & parametres() {return _parametres;}

  //=====================================================//
  //Crop de l'image en prenant des marges
  // imgf         : "grande" image
  // canal        : numero du canal a traiter
  // maxim        : reetaler entre 0 et maxim (si maxim<=0, calcul du max sur la dalle)
  // i0_grande_image, j0_grande_image   : numero ligne et colonne de debut du crop
  // nc, nl       : nombre de lignes et de colonnes maximales souhaitees
  // marge        : marge a prendre en compte (cf pyramide)
  // xmin,ymin,xmax,ymax : Vraie emprise de siftage au sein de la dalle : on ne conservera donc que les points Pt de Lowe  tels que xmin<=Pt.x<=xmax && ymin<=Pt.y<=ymax
  //-----------------------------------------------------//

  const bool GetDalle(const AccesImage & imgf, 
		      const int i0_grande_image, const int j0_grande_image, const int nc, const int nl, 
		      const int marge, int & xmin, int & ymin, int & xmax, int & ymax,
		      VL::ImageBuffer & buffer, 
		      const bool verbose=true) const;

  const bool GetDalle(const AccesImage & imgf, 
		      const ParametragePretraitement parametragepretraitement, 
		      const int i0_grande_image, const int j0_grande_image, const int nc, const int nl, 
		      const int marge, int & xmin, int & ymin, int & xmax, int & ymax,
		      VL::ImageBuffer & buffer, 
		      const bool verbose=true) const;

  const bool GetDalle(const AccesImage & imgf, 
		      const int canal, 
		      const float maxim, 
		      const bool usegradient,
		      const int i0_grande_image, const int j0_grande_image, const int nc, const int nl, 
		      const int marge, int & xmin, int & ymin, int & xmax, int & ymax,
		      VL::ImageBuffer & buffer, 
		      const bool verbose=true) const;
  
  
  //Fonctions appel
  const bool GetDalle(const std::string & adresseimage, 
		      const int canal, 
		      const float maxim,
		      const bool usegradient,
		      const int i0_grande_image, const int j0_grande_image, const int nc, const int nl, 
		      const int marge, int & xmin, int & ymin, int & xmax, int & ymax, 
		      VL::ImageBuffer & buffer, 
		      const bool verbose=true) const;
  
  const bool GetDalle(const std::string & adresseimage, 
		      const ParametragePretraitement parametragepretraitement,
		      const int i0_grande_image, const int j0_grande_image, const int nc, const int nl, 
		      const int marge, int & xmin, int & ymin, int & xmax, int & ymax, 
		      VL::ImageBuffer & buffer, 
		      const bool verbose=true) const;

  /*
  void ChangementFormatExterneVersInterne(const ImageMemoire & im, VL::ImageBuffer & buffer, const int canal=0, const bool usegradient=false, const float maxim=-1, const float minim=0) const;
  void ChangementFormatExterneVersInterne(const ImageMemoire & im, VL::ImageBuffer & buffer, const ParametragePretraitement & parametrage) const;
  void ChangementFormatExterneVersInterne(const ImageMemoire & im, VL::ImageBuffer & buffer) const{
    ChangementFormatExterneVersInterne(im,buffer,_parametres);
  }
  */

  /*
    cIm2D<tREAL4> 
   */

  void ChangementFormatExterneVersInterne(ImageMemoire & im, VL::ImageBuffer & buffer, const int canal=0, const bool usegradient=false, const float maxim=-1, const float minim=0) const;
  void ChangementFormatExterneVersInterne(ImageMemoire & im, VL::ImageBuffer & buffer, const ParametragePretraitement & parametrage) const;
  void ChangementFormatExterneVersInterne(ImageMemoire & im, VL::ImageBuffer & buffer) const{
    ChangementFormatExterneVersInterne(im,buffer,_parametres);
  }

};


#endif
