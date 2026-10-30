#ifndef APPARIEMENT_DESCRIPTEURS_H
#define APPARIEMENT_DESCRIPTEURS_H

#include "KDTree.hpp"


//=============================================================
//Classe AppariementDescripteurs
//=============================================================
class AppariementDescripteurs{
 protected:
  int _maxPtsAvisiter;
  double _ratiodistances;
  bool _normalisedescripteur;
  bool _verbose;

 public:
  //Constructeur / Destructeur
  AppariementDescripteurs(const int maxPtsAvisiter=200, const double ratiodistances=0.6, const bool normalisedescripteur=false, const bool verbose=true):_maxPtsAvisiter(maxPtsAvisiter),_ratiodistances(ratiodistances),_normalisedescripteur(normalisedescripteur),_verbose(verbose){}
  
  //Accesseurs
  const int maxPtsAvisiter() const{return _maxPtsAvisiter;}
  int & maxPtsAvisiter(){return _maxPtsAvisiter;}
  const double ratiodistances() const{return _ratiodistances;}
  double & ratiodistances(){return _ratiodistances;}
  const bool normalisedescripteur() const{return _normalisedescripteur;}
  bool & normalisedescripteur(){return _normalisedescripteur;}
  const bool verbose() const{return _verbose;}
  bool & verbose(){return _verbose;}

  //Fonctions

  //Appariement par la methode SIFT classique
  inline int apparie(const std::string & adresse_kptsdata, const std::string & adresse_kptsquery, const std::string & adresse_result) const;
  inline int apparie(const std::string & adresse_kptsdata, const std::string & adresse_kptsquery, std::vector<T2D<double> > & pts1, std::vector<T2D<double> > & pts2) const;
  

  void apparie(const std::vector<KeyPoint> & kptsdata, const std::vector<KeyPoint> & kptsquery, std::vector<unsigned int> & iddata_app, std::vector<unsigned int> & idquery_app){
    AppariementDescripteurs::Apparie(kptsdata,kptsquery,iddata_app,idquery_app,_maxPtsAvisiter,_ratiodistances,_normalisedescripteur,_verbose);
  }

  void ppv(const std::vector<KeyPoint> & kptsdata, const std::vector<KeyPoint> & kptsquery, std::vector<unsigned int> & iddata_app, std::vector<unsigned int> & idquery_app, std::vector<double> & dist1, std::vector<double> & dist2){
    AppariementDescripteurs::PPV(kptsdata,kptsquery,iddata_app,idquery_app,dist1,dist2,_maxPtsAvisiter,_normalisedescripteur,_verbose);
  }

  //Recherche des nbvois ppv
  void nppv(const std::vector<KeyPoint> & kptsdata, const std::vector<KeyPoint> & kptsquery, const unsigned int nbvois, std::vector<std::vector<unsigned int> > & idptsppv, std::vector<std::vector<double> > & distsppv){
    AppariementDescripteurs::NPPV(kptsdata,kptsquery,nbvois,idptsppv,distsppv,_maxPtsAvisiter,_normalisedescripteur,_verbose);
  }

  //Recherche des nbvois plus proches voisins dans l'espace des descripteurs puis  test avec predicteur et selection parmi eux du plus proche voisin dans l'image s'il se trouve a une distance < tolerance_distance_image
  template<class TPredicteur2D> inline void apparie_nppv_et_predicteur(const std::vector<KeyPoint> & kptsdata, const std::vector<KeyPoint> & kptsquery, const unsigned int nbvois, const TPredicteur2D & predicteur, const double tolerance_distance_image, std::vector<unsigned int> & iddata_app, std::vector<unsigned int> & idquery_app);



  inline static void Apparie(const std::vector<KeyPoint> & kptsdata, const std::vector<KeyPoint> & kptsquery, std::vector<unsigned int> & iddata_app, std::vector<unsigned int> & idquery_app, const int maxPtsAvisiter=200, const double ratiodistances=0.6, const bool normalisedescripteur=false, const bool verbose=true);

  //recherche juste le PPV et renvoie les distances aux 2 PPV...
  inline static void PPV(const std::vector<KeyPoint> & kptsdata, const std::vector<KeyPoint> & kptsquery, std::vector<unsigned int> & iddata_app, std::vector<unsigned int> & idquery_app, std::vector<double> & dist1, std::vector<double> & dist2, const int maxPtsAvisiter=200, const bool normalisedescripteur=false, const bool verbose=true);

//Recherche des nbvois ppv
  inline static void NPPV(const std::vector<KeyPoint> & kptsdata, const std::vector<KeyPoint> & kptsquery, const unsigned int nbvois, std::vector<std::vector<unsigned int> > & idptsppv, std::vector<std::vector<double> > & distsppv, const int maxPtsAvisiter=200, const bool normalisedescripteur=false, const bool verbose=true);

};


void AppariementDescripteurs::Apparie(const std::vector<KeyPoint> & kptsdata, const std::vector<KeyPoint> & kptsquery, std::vector<unsigned int> & iddata_app, std::vector<unsigned int> & idquery_app, const int maxPtsAvisiter, const double ratiodistances, const bool normalisedescripteur, const bool verbose){
  
  //----------------------------------------------------------------------
  //	Parameters
  //----------------------------------------------------------------------
  int				dim				= 128;			// dimension des descripteurs
  double ratiodistancescarrees=pow(ratiodistances,2.0);//pour comparer des distances au carre...


  //Simples verifications avant d'aller plus loin pour rien...
  if(kptsdata.empty()){return;}
  if(kptsquery.empty()){return;}
  std::vector<KeyPoint>::const_iterator ikptsdata=kptsdata.begin();
  std::vector<KeyPoint>::const_iterator ikptsquery=kptsquery.begin();
  if((*ikptsdata).dim()!=(*ikptsquery).dim()){
    std::cout<<"Dimension incompatible entre les points data et query"<<std::endl;
    return;
  }else{dim=(*ikptsdata).dim();}

  iddata_app.clear();
  idquery_app.clear();


  //-------------------------------------------
  //Traitements data et construction du k-d-tree
  //-------------------------------------------
  //Allocations...

  KDTree kdtree;
  kdtree.set_size_dim(kptsdata.size(),dim);
  

  //Lecture des points data
  if(verbose)std::cout<<"Chargement des points data..."<<std::flush;
  kdtree.set_data(kptsdata);
  if(verbose)std::cout<<"OK : "<<kptsdata.size()<<" points data"<<std::endl;

  //Construction du kd-tree
  if(verbose)std::cout<<"Construction du k-d-tree..."<<std::flush;
  kdtree.construire_kdtree();
  if(verbose)std::cout<<"OK"<<std::endl;


  //-------------------------------------------
  //Lecture des points query et recherche de leurs NN dans l'arbre
  //-------------------------------------------
  if(verbose)std::cout<<"Chargement des points query et requetes..."<<std::flush;
  if(verbose)std::cout<<kptsquery.size()<<" requetes a traiter..."<<std::flush;

  //Allocation
  std::vector<unsigned int> nnIdx;
  std::vector<double> dists;

  ikptsquery=kptsquery.begin();
    int iquery=0;
    for(;ikptsquery!=kptsquery.end();++ikptsquery,iquery++){
    kdtree.nppv((*ikptsquery),2,nnIdx,dists,maxPtsAvisiter,false);

    if(dists[0]<ratiodistancescarrees*dists[1]){
      iddata_app.push_back(nnIdx[0]);
      idquery_app.push_back(iquery);
    }
  }
  if(verbose)std::cout<<"OK"<<std::endl;
  
}


void AppariementDescripteurs::PPV(const std::vector<KeyPoint> & kptsdata, const std::vector<KeyPoint> & kptsquery, std::vector<unsigned int> & iddata_app, std::vector<unsigned int> & idquery_app, std::vector<double> & dist1, std::vector<double> & dist2, const int maxPtsAvisiter, const bool normalisedescripteur, const bool verbose){
  
  //----------------------------------------------------------------------
  //	Parameters
  //----------------------------------------------------------------------
  int				dim				= 128;			// dimension des descripteurs


  //Simples verifications avant d'aller plus loin pour rien...
  if(kptsdata.empty()){return;}
  if(kptsquery.empty()){return;}
  std::vector<KeyPoint>::const_iterator ikptsdata=kptsdata.begin();
  std::vector<KeyPoint>::const_iterator ikptsquery=kptsquery.begin();
  if((*ikptsdata).dim()!=(*ikptsquery).dim()){
    std::cout<<"Dimension incompatible entre les points data et query"<<std::endl;
    return;
  }else{dim=(*ikptsdata).dim();}

  iddata_app.clear();
  idquery_app.clear();
  dist1.clear();
  dist2.clear();

  //-------------------------------------------
  //Traitements data et construction du k-d-tree
  //-------------------------------------------
  //Allocations...

  KDTree kdtree;
  kdtree.set_size_dim(kptsdata.size(),dim);
  

  //Lecture des points data
  if(verbose)std::cout<<"Chargement des points data..."<<std::flush;
  kdtree.set_data(kptsdata);
  if(verbose)std::cout<<"OK : "<<kptsdata.size()<<" points data"<<std::endl;

  //Construction du kd-tree
  if(verbose)std::cout<<"Construction du k-d-tree..."<<std::flush;
  kdtree.construire_kdtree();
  if(verbose)std::cout<<"OK"<<std::endl;


  //-------------------------------------------
  //Lecture des points query et recherche de leurs NN dans l'arbre
  //-------------------------------------------
  if(verbose)std::cout<<"Chargement des points query et requetes..."<<std::flush;
  if(verbose)std::cout<<kptsquery.size()<<" requetes a traiter..."<<std::flush;

  //Allocation
  std::vector<unsigned int> nnIdx;
  std::vector<double> dists;

  ikptsquery=kptsquery.begin();
  unsigned int iquery=0;
  for(;ikptsquery!=kptsquery.end();iquery++,++ikptsquery){
    if(normalisedescripteur)(*ikptsquery).normalise();
    kdtree.nppv((*ikptsquery),2,nnIdx,dists,maxPtsAvisiter,false);
    
    iddata_app.push_back(nnIdx[0]);
    idquery_app.push_back(iquery);
    dist1.push_back(dists[0]);
    dist2.push_back(dists[1]);
  }
  if(verbose)std::cout<<"OK"<<std::endl;
}

void AppariementDescripteurs::NPPV(const std::vector<KeyPoint> & kptsdata, const std::vector<KeyPoint> & kptsquery, const unsigned int nbppv, std::vector<std::vector<unsigned int> > & idptsppv, std::vector<std::vector<double> > & distsppv, const int maxPtsAvisiter, const bool normalisedescripteur, const bool verbose){
  
  //----------------------------------------------------------------------
  //	Parameters
  //----------------------------------------------------------------------
  int				dim				= 128;			// dimension des descripteurs


  //Simples verifications avant d'aller plus loin pour rien...
  if(kptsdata.empty()){return;}
  if(kptsquery.empty()){return;}
  std::vector<KeyPoint>::const_iterator ikptsdata=kptsdata.begin();
  std::vector<KeyPoint>::const_iterator ikptsquery=kptsquery.begin();
  if((*ikptsdata).dim()!=(*ikptsquery).dim()){
    std::cout<<"Dimension incompatible entre les points data et query"<<std::endl;
    return;
  }else{dim=(*ikptsdata).dim();}

  idptsppv.clear();
  distsppv.clear();

  idptsppv.resize(kptsquery.size(),std::vector<unsigned int>(nbppv));
  distsppv.resize(kptsquery.size(),std::vector<double>(nbppv));

  //-------------------------------------------
  //Traitements data et construction du k-d-tree
  //-------------------------------------------
  //Allocations...

  KDTree kdtree;
  kdtree.set_size_dim(kptsdata.size(),dim);
  

  //Lecture des points data
  if(verbose)std::cout<<"Chargement des points data..."<<std::flush;
  kdtree.set_data(kptsdata);
  if(verbose)std::cout<<"OK : "<<kptsdata.size()<<" points data"<<std::endl;

  //Construction du kd-tree
  if(verbose)std::cout<<"Construction du k-d-tree..."<<std::flush;
  kdtree.construire_kdtree();
  if(verbose)std::cout<<"OK"<<std::endl;


  //-------------------------------------------
  //Lecture des points query et recherche de leurs NN dans l'arbre
  //-------------------------------------------
  if(verbose)std::cout<<"Chargement des points query et requetes..."<<std::flush;
  if(verbose)std::cout<<kptsquery.size()<<" requetes a traiter"<<std::endl;

  std::vector<unsigned int> nnIdx;
  std::vector<double> dists;

  std::vector<std::vector<unsigned int> >::iterator i_idpts=idptsppv.begin();
  std::vector<std::vector<double> >::iterator i_dist=distsppv.begin();

  ikptsquery=kptsquery.begin();
  unsigned int iquery=0;
  for(;ikptsquery!=kptsquery.end();iquery++,++ikptsquery,++i_idpts,++i_dist){
    kdtree.nppv((*ikptsquery),nbppv,nnIdx,dists,maxPtsAvisiter,false);
    (*i_idpts)=nnIdx;
    (*i_dist)=dists;    
  }
  if(verbose)std::cout<<"OK"<<std::endl;

}


//----------------------------------------------------------------------
//Recherche des nbvois plus proches voisins dans l'espace des descripteurs puis test avec predicteur et selection parmi eux du plus proche voisin dans l'image s'il se trouve a une distance < tolerance_distance_image
//----------------------------------------------------------------------
template<class TPredicteur2D>
inline void AppariementDescripteurs::apparie_nppv_et_predicteur(const std::vector<KeyPoint> & kptsdata, const std::vector<KeyPoint> & kptsquery, const unsigned int nbvois, const TPredicteur2D & predicteur, const double tolerance_distance_image, std::vector<unsigned int> & iddata_app, std::vector<unsigned int> & idquery_app){
  double tolerance_distance2_image=std::pow(tolerance_distance_image,2.0);//On travaille avec des distances au carre...


  //-------------------------------------------
  //Appariement multiple : recherche des nbvois plus proches voisins descripteur
  //-------------------------------------------
  std::vector<std::vector<unsigned int> > idpts;
  std::vector<std::vector<double> > dist;
  if(_verbose){std::cout<<"Recherche des plus proches voisins descripteur..."<<std::endl;}
  nppv(kptsdata,kptsquery,nbvois,idpts,dist);
    
  //-------------------------------------------
  //Test croise avec predicteur : on ne conserve que le plus proche voisin image s'il se trouve a une distance < tolerance_distance_image
  //-------------------------------------------
  if(_verbose){std::cout<<"Recherche des plus proches voisins image..."<<std::flush;}
  iddata_app.reserve(idpts.size());
  idquery_app.reserve(idpts.size());
  std::vector<std::vector<unsigned int> >::const_iterator i_idpts=idpts.begin();
  T2D<double> ptpred,pttmp;
  double disttmp,dmin;
  int kmin;
  for(int n=0;n<idpts.size();n++,++i_idpts){
    ptpred=T2D<double>(predicteur(kptsquery[n].xy()));
    pttmp=T2D<double>(kptsdata[(*i_idpts)[0]].xy());
    disttmp=std::pow(ptpred.x()-pttmp.x(),2.0)+std::pow(ptpred.y()-pttmp.y(),2.0);
    dmin=disttmp;
    kmin=0;
    for(int k=1;k<nbvois;k++){
      pttmp=T2D<double>(kptsdata[(*i_idpts)[k]].xy());
      disttmp=std::pow(ptpred.x()-pttmp.x(),2.0)+std::pow(ptpred.y()-pttmp.y(),2.0);
      if(disttmp<dmin){dmin=disttmp; kmin=k;}
    }
    if(dmin>tolerance_distance2_image){continue;}
    iddata_app.push_back((*i_idpts)[kmin]);
    idquery_app.push_back(n);
  }
  if(_verbose){std::cout<<"OK"<<std::endl;}
  
}



#endif
