#ifndef _SIFTDRIVER_HPP_
#define _SIFTDRIVER_HPP_

#include "sift.hpp"
#include "KeyPoint.h"
#include "AccesImage.hpp"
#include "PreparationImage.hpp"
#include "sift-driver_parametres.hpp"
#include "outilsstring/Stringpp.hpp"


//=========================================================
//SiftDriver
//Classe de surcouche par rapport aux classes de sift++ (Vedaldi)
// 
// Le parametrage est assure par ParametrageSiftDriver
// Permet d'extaire les points de Lowe d'une image en la dallant
//=========================================================
class SiftDriver{
 protected:
  ParametrageSiftDriver _parametres;

 public:
  //------------------------------
  //Constructeurs - Destructeurs
  SiftDriver(const ParametrageSiftDriver & parametres):_parametres(parametres){}
  SiftDriver(){}
  virtual ~SiftDriver(){}

  //------------------------------
  //Acces aux parametres
  const ParametrageSiftDriver & parametres() const{return _parametres;}
  ParametrageSiftDriver & parametres() {return _parametres;}


  //------------------------------
  //Fonctions de traitement

  //Traitement d'une image complete
  template<class TFluxSortieKeypoint>
    inline const int traiter_image(const std::string & adresseimg,
				   TFluxSortieKeypoint & ffkey);
  
  //Traitement de parties { [ coordimage_coincrop[n].x();coordimage_coincrop[n].x()+taillecrop[n].x() [ x [ coordimage_coincrop[n].y();coordimage_coincrop[n].y()+taillecrop[n].y() [ } d'une image
  template<class TFluxSortieKeypoint>
    inline const int traiter_image(const std::string & adresseimg,
				   const std::vector<T2D<double> > & coordimage_coincrop,
				   const std::vector<T2D<double> > & taillecrop,
				   TFluxSortieKeypoint & ffkey);

  //Traitement d'une partie [X0;X0+NbColonnes[ x [Y0;Y0+NbLignes[ d'une image
  template<class TFluxSortieKeypoint>
    inline const int traiter_image(const std::string & adresseimg,
				   const double X0, const double Y0,
				   const double NbColonnes, const double NbLignes,
				   TFluxSortieKeypoint & ffkey);
  
  
  //Traitement d'une image complete
  template<class TFluxSortieKeypoint>
    inline const int traiter_image(const AccesImage & imgf,
				   TFluxSortieKeypoint & ffkey);
  
  //Traitement de parties { [ coordimage_coincrop[n].x();coordimage_coincrop[n].x()+taillecrop[n].x() [ x [ coordimage_coincrop[n].y();coordimage_coincrop[n].y()+taillecrop[n].y() [ } d'une image
  template<class TFluxSortieKeypoint>
    inline const int traiter_image(const AccesImage & imgf,
				   const std::vector<T2D<double> > & coordimage_coincrop,
				   const std::vector<T2D<double> > & taillecrop,
				   TFluxSortieKeypoint & ffkey);
  
  //Traitement d'une partie [X0;X0+NbColonnes[ x [Y0;Y0+NbLignes[ d'une image
  template<class TFluxSortieKeypoint>
    inline const int traiter_image(const AccesImage & imgf,
				   const double X0, const double Y0,
				   const double NbColonnes, const double NbLignes,
				   TFluxSortieKeypoint & ffkey);


  //--------------------------------
  //Operation elementaire :Traitement d'une dalle
  // buffer : image (dans le format de Vedaldi)
  // xmin,ymin,xmax,ymax : emprise "utile" au sein de la dalle (notamment si l'on utilise des marges) ; on ne conservera que les points qui s'y trouve
  // X0,Y0 : coordonnees du coin superieur gauche de la dalle dans la grande image
  // fkey : flux de sortie
  template<class TFluxSortieKeypoint> 
    inline const int traiter_dalle(VL::ImageBuffer & buffer,
			    const int xmin, const int ymin, const int xmax, const int ymax,//emprise "utile"
			    const double X0, const double Y0,//coordonnees coin superieur gauche
			    TFluxSortieKeypoint & fkey);
  

};



template<class TFluxSortieKeypoint>
inline const int SiftDriver::traiter_image(const std::string & adresseimg,
					   TFluxSortieKeypoint & ffkey){
  return traiter_image(AccesImage(adresseimg),ffkey);
}

template<class TFluxSortieKeypoint>
inline const int SiftDriver::traiter_image(const std::string & adresseimg,
					   const std::vector<T2D<double> > & coordimage_coincrop,
					   const std::vector<T2D<double> > & taillecrop,
					   TFluxSortieKeypoint & ffkey){
  return traiter_image(AccesImage(adresseimg),coordimage_coincrop,taillecrop,ffkey);
}
 

template<class TFluxSortieKeypoint>
inline const int SiftDriver::traiter_image(const std::string & adresseimg,
					   const double X0, const double Y0,
					   const double NbColonnes, const double NbLignes,
					   TFluxSortieKeypoint & ffkey){
  return traiter_image(AccesImage(adresseimg),X0,Y0,NbColonnes,NbLignes,ffkey);
}
 
 
template<class TFluxSortieKeypoint>
inline const int SiftDriver::traiter_image(const AccesImage & imgf,
					   TFluxSortieKeypoint & ffkey){
  return traiter_image(imgf,0,0,imgf.GetTaille().X(),imgf.GetTaille().Y(),ffkey);
}

 
template<class TFluxSortieKeypoint>
inline const int SiftDriver::traiter_image(const AccesImage & imgf,
					   const std::vector<T2D<double> > & coordimage_coincrop,
					   const std::vector<T2D<double> > & taillecrop,
					   TFluxSortieKeypoint & ffkey){
  if(coordimage_coincrop.size()!=taillecrop.size()){
    std::cout<<"Pb : pas meme nb crop..."<<std::endl;
    return 0;
  }
  std::vector<T2D<double> >::const_iterator it_coincrop=coordimage_coincrop.begin();
  std::vector<T2D<double> >::const_iterator it_taillecrop=taillecrop.begin();
  for(;it_coincrop!=coordimage_coincrop.end();++it_coincrop,++it_taillecrop){
     traiter_image(imgf,(*it_coincrop).x(),(*it_coincrop).y(),(*it_taillecrop).x(),(*it_taillecrop).y(),ffkey);
  }
}


template<class TFluxSortieKeypoint>
inline const int SiftDriver::traiter_image(const AccesImage & imgf,
					   const double X0, const double Y0,
					   const double NbColonnes, const double NbLignes,
					   TFluxSortieKeypoint & ffkey){
  
  //On commence par verifier si les parametres fournis sont coherents avec l'image...
  double Nc=imgf.GetTaille().X();
  double Nl=imgf.GetTaille().Y();
  if(X0>=Nc || Y0>=Nl){return 0;}
  if(NbColonnes<=0 || NbLignes<=0){return 0;}
  double _X0=X0;
  double _Y0=Y0;
  double _NbColonnes=NbColonnes;
  double _NbLignes=NbLignes;
  if(_X0<0){
    _X0=0.0;
    _NbColonnes=_NbColonnes-_X0;
  }      
  if(_Y0<0){
    _Y0=0.0;
    _NbLignes=_NbLignes-_Y0;
  }      
  _NbColonnes=std::min(Nc,_NbColonnes);
  _NbLignes=std::min(Nl,_NbLignes);
  
  //Calcul du decoupage en dalles
  Dallage dallage(parametres().dallage());
  dallage.verbose()=parametres().verbose();
  
  std::vector<int> nc_dalles,nl_dalles;
  dallage.ChoixMeilleursParametresDecoupe(_NbColonnes,nc_dalles);//en colonnes
  dallage.ChoixMeilleursParametresDecoupe(_NbLignes,nl_dalles);//en lignes
  
  //Boucle sur les images decoupees
  T2D<double> coordimage_coincrop(_X0,_Y0);
  for(std::vector<int>::const_iterator ic=nc_dalles.begin();ic!=nc_dalles.end();++ic){
    coordimage_coincrop.y()=_Y0;
    for(std::vector<int>::const_iterator il=nl_dalles.begin();il!=nl_dalles.end();++il){
      VL::ImageBuffer buffer ;
      int xmin, ymin, xmax, ymax;//emprise utile (i.e. hors marge) dans le crop
      //std::cout<<coordimage_coincrop.x()<<" "<<coordimage_coincrop.y()<<std::endl;continue;
      if(!PreparationImage().GetDalle(imgf,
				      parametres().pretraitement(),
				      coordimage_coincrop.x(),coordimage_coincrop.y(),(*ic),(*il),parametres().dallage().marge(),
				      xmin,ymin,xmax,ymax,
				      buffer,parametres().verbose())
	 ){continue;}
      /*int reussite=*/traiter_dalle(buffer,
				 xmin,ymin,xmax,ymax,
				 coordimage_coincrop.x(),coordimage_coincrop.y(),
				 ffkey);
      coordimage_coincrop.y()=coordimage_coincrop.y()+(*il);
    }
    coordimage_coincrop.x()=coordimage_coincrop.x()+(*ic);
  }
  return 0;
}
 

template<class TFluxSortieKeypoint> 
inline const int SiftDriver::traiter_dalle(VL::ImageBuffer & buffer,
					   const int xmin, const int ymin, const int xmax, const int ymax,//emprise "utile"
					   const double X0, const double Y0,//coordonnees coin superieur gauche
					   TFluxSortieKeypoint & fkey){
  
  int omin=parametres().sift().firstoctave();
  int O=parametres().sift().nboctaves();
  int S=parametres().sift().nblevels();
  
  // ---------------------------------------------------------------
  //                                            Gaussian scale space
  // ---------------------------------------------------------------    
  if(parametres().verbose())std::cout << "Calcul Gaussian scale space ... " <<std::flush ;
  
  // optionally autoselect the number number of octaves
  if(O < 1) {
    O=std::max(int(std::floor(log2(std::min(buffer.width,buffer.height)))-omin-4), 1) ;
  }
  
  
  // initialize scalespace
  VL::Sift sift(buffer.data, buffer.width, buffer.height, 
		parametres().sift().sigman(), parametres().sift().sigma0(),
		O, S,
		omin, -1, S+1) ;
  
  parametres().verbose() && std::cout << "("<<O<<" octaves a partir de "<<omin<<", "<<S<< " niveaux par octave) OK."<<std::endl ;
  
  delete [] buffer.data;
  
  /*
  // ---------------------------------------------------------------
  //                                       Save Gaussian scale space
  // ---------------------------------------------------------------    
  
  if(savegss) {
  verbose &&std:: cout<<"Saving Gaussian scale space."<<std::endl ;
    stringpp imageBasename(outputFilename.newSupprExtension());
    for(int o = omin ; o < omin + O ; ++o) {
      for(int s = 0 ; s < S ; ++s) {
        std::ostringstream suffix ;
        suffix<<'.'<<o<<'.'<<s<<".pgm" ;
        std::string imageFilename = imageBasename + suffix.str() ;
        verbose && std::cout << "Saving octave " <<o
                             << " level "<<s
	       		     << " to '"<<imageFilename<<"' ..."<<std::flush ;
        std::ofstream fout(imageFilename.c_str()) ;
        if(!fout.good()) 
	      VL_THROW("Could not open '"<<imageFilename<<'\'') ;
	      
        VL::insertPgm(fout,sift.getLevel(o,s),sift.getOctaveWidth(o),sift.getOctaveHeight(o));
        fout.close();
        verbose && std::cout<<" done."<<std::endl ;
      }
    }
  }
  */
      
  // ---------------------------------------------------------------
  //                                                   SIFT detector
  // ---------------------------------------------------------------    
  
  if(parametres().sift().haveKeypoints()==0 && parametres().sift().siftdense()==0) {
    parametres().verbose() && std::cout << "Detecteur SIFT (threshold:"<<parametres().sift().threshold()
					<< ", edge-threshold: "<<parametres().sift().edgeThreshold()
					<< ") ..."<<std::flush ;
    
    sift.detectKeypoints(parametres().sift().threshold(), parametres().sift().edgeThreshold(), parametres().sift().signe(), xmin,ymin,xmax,ymax);

    parametres().verbose() && std::cout <<"OK "<<std::endl;
  }else if(parametres().sift().haveKeypoints()!=0){
    //Reads keypoint from file, compute descriptors
    //A REACTIVER UN JOUR, mais plutot plus haut (cf lecture et decoupage en dalle)
    for(std::vector<T3D<float> >::const_iterator ipadecrire=_parametres.sift().KeyPointsADecrire().begin();ipadecrire!=_parametres.sift().KeyPointsADecrire().end();++ipadecrire){
      T3D<float> kptmp((*ipadecrire).X(),(*ipadecrire).Y(),(*ipadecrire).Z());
      kptmp.x()=kptmp.X()-X0+xmin;
      kptmp.y()=kptmp.Y()-X0+ymin;
      if(kptmp.x()<xmin || kptmp.y()<ymin || kptmp.x()>xmax || kptmp.y()>ymax){continue;}
      sift.setKeypoints(sift.getKeypoint(kptmp.X(),kptmp.Y(),kptmp.Z()));
    }
  }else if(parametres().sift().issiftdense()){//if(parametres().sift().siftdense()!=0){
    double iz;
    if(parametres().sift().issiftdensescalefixe()){
      iz=parametres().sift().siftdensescalefixe(); 
    }else{
      iz=sift.getScaleFromIndex(parametres().sift().firstoctave(),0);
    }
    //std::cout<<iz<<std::endl;
    //if(verbose){std::cout<<iz<<std::endl;}
    //double iz=1.0;//ca marche pour les vignes...
    for(double ix=xmin;ix<xmax;ix=ix+parametres().sift().siftdense()){
      for(double iy=ymin;iy<ymax;iy=iy+parametres().sift().siftdense()){
	sift.setKeypoints(sift.getKeypoint(ix,iy,iz));
      } 
    }
  }
  
  // ---------------------------------------------------------------
  //                               SIFT orientations and descriptors
  // ---------------------------------------------------------------    
  
  if( parametres().verbose() ) {
    if( !parametres().sift().noorient() &&  parametres().sift().nodescr()) std::cout<<"Calcul orientations...";
    if(  parametres().sift().noorient() && !parametres().sift().nodescr()) std::cout<<"Calcul descripteurs...";
    if( !parametres().sift().noorient() && !parametres().sift().nodescr()) std::cout<<"Calcul orientations and descripteurs...";
    if(  parametres().sift().noorient() &&  parametres().sift().nodescr()) std::cout<<"Fin..."; 
    std::cout<<std::flush ;
  }
  
  {     
    int taille_descr=128;
    if(parametres().sift().nodescr()){taille_descr=0;}
    KeyPoint kp(taille_descr);
    VL::float_t descr[128];
    {
      // -------------------------------------------------------------
      // Compute orientations and descriptors
      // -------------------------------------------------------------
      VL::float_t angles[4] ;
      for(VL::Sift::KeypointsConstIter iter=sift.keypointsBegin();iter!=sift.keypointsEnd();++iter){
	//std::cout<<iter->o<<" "<<iter->x<<" "<<iter->y<<" "<<iter->s<<" "<<iter->ix<<" "<<iter->iy<<" "<<iter->is<<" "<<iter->sigma<<std::endl;//affichage test 29/06/2011

	//Prise en compte des marges (et a garder dans le cas general)
	if(iter->x<xmin || iter->y<ymin || iter->x>xmax || iter->y>ymax){
	  //std::cout<<"HZ : "<<iter->x<<" "<<iter->y<<std::endl;
	   continue;
	}
	
	// detect orientations
	int nangles ;
	if( ! parametres().sift().noorient() ) {
	  nangles = sift.computeKeypointOrientations(angles, *iter) ;
	}else{
	  nangles = 1;
	  angles[0] = VL::float_t(0) ;
	}
	
	// compute descriptors
	for(int a = 0 ; a < nangles ; ++a) {
	  kp.x()=iter->x+X0-xmin;
	  kp.y()=iter->y+Y0-ymin;
	  kp.sigma()=iter->sigma;
	  kp.angle()=angles[a];
	  kp.signe()=iter->signe;//Ajout du "signe" du point (< ou >)
	  
	  if(!parametres().sift().nodescr()) {
	    sift.computeKeypointDescriptor(descr, *iter, angles[a]);
	    KeyPoint::iterator idesc=kp.begin();
	    for(int i=0;i<128;++i,++idesc){(*idesc)=VL::uint8_t(VL::float_t(512)*descr[i]);}
	    if(parametres().sift().demidescripteur()){kp.DemiDescripteur(8,true);}
	  }
	  fkey<<kp;
	  if(parametres().sift().demidescripteur()){kp.desc().resize(128);}//sinon pb pour le point suivant...
	  
	} // next angle
      } // next keypoint
    }
    //continue;//modif test memoire
    
    parametres().verbose() && std::cout << " OK."<<std::endl ;
  }
 
  return 0 ;
}


#endif
