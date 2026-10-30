#ifndef _KDTREE_HPP_
#define _KDTREE_HPP_


//================================================================================
// Classe KDTree
//      classe encapsulant le kdtree de ANN
//      permet de mieux gérer le chargement des points (en une fois, ou en plusieurs...)
//      constructeur, destructeur --> gestion des allocations memoire
//--------------------------------------------------------------------------------

class KDTree{

protected:
  //Allocations...
  bool _data_ok;
  int _i_data;
  bool _kdtree_ok;
  int _dim;
  int			_ndataPts;				// actual number of data points
  ANNpointArray		_dataPts;				// data points
  ANNkd_tree*		_kdTree;				// search structure
  bool _normalisedescripteur;

public:
  KDTree():_data_ok(false),_i_data(0),_kdtree_ok(false),_dim(0),_ndataPts(0),_normalisedescripteur(false){}

  virtual ~KDTree(){
    if(_kdtree_ok){
      delete _kdTree;
    }
    if(_data_ok){
      annDeallocPts(_dataPts);
    }
    annClose();									// done with ANN
  }

  //----------------------------------------------------------
  //Initialisation du tableau de donnees
  //----------------------------------------------------------
  void set_size_dim(const int size, const int dim){
    if(_data_ok){return;}//deja fait, on ne doit plus modifier
    _ndataPts=size;
    _dim=dim;
    if(_dim==0 || _ndataPts==0){return;}
    _dataPts=annAllocPts(_ndataPts,_dim);			// allocate data points
    _data_ok=true;
    _i_data=0;
    //std::cout<<_ndataPts<<" "<<_dim<<std::endl; return;
  }
  void set_dim(const int dim){set_size_dim(_ndataPts,dim);}
  void set_size(const int size){set_size_dim(size,_dim);}

  //Simples accesseurs
  const unsigned int dim() const {return _dim;}
  const unsigned int size() const {return _ndataPts;}
  ANNkd_tree* kdTree(){return _kdTree;}


  //----------------------------------------------------------
  //Chargement des donnees
  //Dans ce cas de figure, on charge toutes les donnees en une fois... 
  //----------------------------------------------------------
  //donnees quelconques
  template<class T>
  void set_data(const std::vector<std::vector<T> > & data){
    if(data.empty()){return;}
    int size=data.size();
    int dim=data.front().size();
    if(dim==0){return;}
    set_size_dim(size,dim);
    typename std::vector<std::vector<T> >::const_iterator itdata=data.begin();
    _i_data=0;
    for(itdata=data.begin();itdata!=data.end();_i_data++,++itdata){
      for(int d=0;d<dim;d++){
	_dataPts[_i_data][d]=(*itdata)[d];
      }
    }
    

  }
  
#ifdef _KEYPOINT_HPP_
  //descripteurs SIFT
  void set_data(const std::vector<KeyPoint> & data){
    if(data.empty()){return;}
    int size=data.size();
    int dim=data.front().dim();
    if(dim==0){return;}
    set_size_dim(size,dim);
    std::vector<KeyPoint>::const_iterator itdata=data.begin();
    _i_data=0;
    for(itdata=data.begin();itdata!=data.end();_i_data++,++itdata){
      for(int d=0;d<dim;d++){
	_dataPts[_i_data][d]=(*itdata)[d];
      }
    }
  }
#endif //_KEYPOINT_HPP_
  
  //----------------------------------------------------------
  //Chargement des donnees
  //Dans ce cas de figure, on charge fichier par fichier
  //...mais il faut tout de même avoir fourni le nombre TOTAL de points auparavant avec set_size_dim() (cf allocation)
  //----------------------------------------------------------
  //donnees quelconques
  template<class T>
  const bool add_data(const std::vector<std::vector<T> > & data){
    int dim=data.front().size();
    if(_dim!=(*this).dim()){return false;}//dimensions incompatibles
    typename std::vector<std::vector<T> >::const_iterator itdata=data.begin();
    for(itdata=data.begin();itdata!=data.end();_i_data++,++itdata){				// read the data
      //if(_normalisedescripteur)(*itdata).normalise();
      for(int d=0;d<dim;d++){
	_dataPts[_i_data][d]=(*itdata)[d];
      }
    }
    return true;
  }
  

#ifdef _KEYPOINT_HPP_
 //descripteurs SIFT
  const bool add_data(const std::vector<KeyPoint> & data){
    int dim=data.front().dim();
    if(_dim!=(int)(*this).dim()){return false;}//dimensions incompatibles
    std::vector<KeyPoint>::const_iterator itdata=data.begin();
    for(itdata=data.begin();itdata!=data.end();_i_data++,++itdata){				// read the data
      //if(_normalisedescripteur)(*itdata).normalise();
      for(int d=0;d<dim;d++){
	_dataPts[_i_data][d]=(*itdata)[d];
      }
    }
    return true;
  }
#endif//_KEYPOINT_HPP_

  //----------------------------------------------------------
  //Chargement des donnees
  //Dans ce cas de figure, on charge fichier par fichier
  //...mais il faut tout de même avoir fourni le nombre TOTAL de points auparavant avec set_size_dim() (cf allocation)
  //----------------------------------------------------------
  //donnees quelconques
  template<class T>
  const bool add_data(const std::vector<T> & data){
    int dim=data.size();
    if(_dim!=(*this).dim()){return false;}//dimensions incompatibles
    for(int d=0;d<dim;d++){
      _dataPts[_i_data][d]=data[d];
    }
    _i_data++;
    return true;
  }
    
#ifdef _KEYPOINT_HPP_
 //descripteurs SIFT
   const bool add_data(const KeyPoint & data){
    int dim=data.dim();
    if(_dim!=(int)(*this).dim()){return false;}//dimensions incompatibles
    for(int d=0;d<dim;d++){
      _dataPts[_i_data][d]=data[d];
    }
    _i_data++;
    return true;
  }
#endif//_KEYPOINT_HPP_


  //----------------------------------------------------------
  //Construction du kdtree a partir des points charges au prealable
  //----------------------------------------------------------
  void construire_kdtree(){
    if(_kdtree_ok){return;}//...deja fait, on ne doit plus modifier
    _kdTree = new ANNkd_tree(				// build search structure
			    _dataPts,			// the data points
			    _ndataPts,			// number of points
			    _dim);			// dimension of space
    _kdtree_ok=true;    
  }


  //----------------------------------------------------------
  //Requetes de n plus proches voisins.
  //----------------------------------------------------------
  //recherche des nppv
#ifdef _KEYPOINT_HPP_
  void nppv(const std::vector<KeyPoint> & kptsquery, const unsigned int nbvois, std::vector<std::vector<unsigned int> > & idptsppv, std::vector<std::vector<double> > & distsppv, const int maxPtsAvisiter=200, const bool verbose=true){
    int				nbnn				= nbvois;			// number of nearest neighbors
    double	       		eps				= 0;			// error bound

    
    if(kptsquery.empty()){return;}
    std::vector<KeyPoint>::const_iterator ikptsquery=kptsquery.begin();
    if(dim()!=(unsigned int)(*ikptsquery).dim()){
      std::cout<<"Dimension incompatible entre les points data et query"<<std::endl;
      return;
    }

    idptsppv.clear();
    distsppv.clear();
    
    idptsppv.resize(kptsquery.size(),std::vector<unsigned int>(nbvois));
    distsppv.resize(kptsquery.size(),std::vector<double>(nbvois));
    

    //-------------------------------------------
    //Lecture des points query et recherche de leurs NN dans l'arbre
    //-------------------------------------------
    if(verbose)std::cout<<"OK\nChargement des points query et requetes..."<<std::flush;
    if(verbose)std::cout<<kptsquery.size()<<" requetes a traiter"<<std::endl;

    //Allocation
    ANNpoint		queryPt;				// query point
    queryPt = annAllocPt(dim());	// allocate query point
    ANNidxArray		nnIdx;					// near neighbor indices
    ANNdistArray	        dists;					// near neighbor distances
    nnIdx = new ANNidx[nbnn];				// allocate near neigh indices
    dists = new ANNdist[nbnn];				// allocate near neighbor dists
    
    std::vector<std::vector<unsigned int> >::iterator i_idpts=idptsppv.begin();
    std::vector<std::vector<double> >::iterator i_dist=distsppv.begin();
    
    ikptsquery=kptsquery.begin();
    unsigned int iquery=0;
    for(;ikptsquery!=kptsquery.end();iquery++,++ikptsquery,++i_idpts,++i_dist){
      //if(normalisedescripteur)(*ikptsquery).normalise();
      for(unsigned int d=0;d<dim();d++){
	queryPt[d]=(*ikptsquery)[d];
      }
    
      annMaxPtsVisit(maxPtsAvisiter);//permet de limiter la recherche aux maxPtsAvisiter (par defaut 200) premiers points visités
      //kdTree->annkSearch(						// search
      _kdTree->annkPriSearch(
			     queryPt,		// query point
			     nbnn,			// number of near neighbors
			     nnIdx,			// nearest neighbors (returned)
			     dists,			// distance (returned)
			     eps);			// error bound
      
      for(int k=0;k<(int)nbvois;k++){
	(*i_idpts)[k]=nnIdx[k];
	(*i_dist)[k]=dists[k];
      }
    }
    if(verbose)std::cout<<"OK"<<std::endl;
    
    delete [] nnIdx;							// clean things up
    delete [] dists;
    annDeallocPt(queryPt);

}



  //recherche des nppv
  void nppv(const KeyPoint & kptquery, const unsigned int nbvois, std::vector<unsigned int> & idptsppv, std::vector<double> & distsppv, const int maxPtsAvisiter, const bool verbose=true){
    nppv(kptquery.desc(),nbvois,idptsppv,distsppv,maxPtsAvisiter,verbose);
  }
#endif //_KEYPOINT_HPP_


  template<class T>
  void nppv(const std::vector<T> & kptquery, const unsigned int nbvois, std::vector<unsigned int> & idptsppv, std::vector<double> & distsppv, const int maxPtsAvisiter, const bool verbose=true){
    int				nbnn				= nbvois;			// number of nearest neighbors
    double	       		eps				= 0;			// error bound

    
    if(dim()!=kptquery.size()){
      std::cout<<"Dimension incompatible entre les points data et query"<<std::endl;
      return;
    }

    
    idptsppv.resize(nbvois);
    distsppv.resize(nbvois);
    

    //-------------------------------------------
    //Lecture des points query et recherche de leurs NN dans l'arbre
    //-------------------------------------------
    if(verbose)std::cout<<"OK\nRequete..."<<std::flush;

    //Allocation
    ANNpoint		queryPt;				// query point
    queryPt = annAllocPt(dim());	// allocate query point
    ANNidxArray		nnIdx;					// near neighbor indices
    ANNdistArray	        dists;					// near neighbor distances
    nnIdx = new ANNidx[nbnn];				// allocate near neigh indices
    dists = new ANNdist[nbnn];				// allocate near neighbor dists
    
    //if(normalisedescripteur)(*kptquery).normalise();
    for(unsigned int d=0;d<dim();d++){
      queryPt[d]=kptquery[d];
    }
    
    annMaxPtsVisit(maxPtsAvisiter);//permet de limiter la recherche aux maxPtsAvisiter (par defaut 200) premiers points visités
    //kdTree->annkSearch(						// search
    _kdTree->annkPriSearch(
			  queryPt,		// query point
			  nbnn,			// number of near neighbors
			  nnIdx,			// nearest neighbors (returned)
			  dists,			// distance (returned)
			  eps);			// error bound
    
    for(unsigned int k=0;k<nbvois;k++){
      idptsppv[k]=nnIdx[k];
      distsppv[k]=dists[k];
    }
    if(verbose)std::cout<<"OK"<<std::endl;
    
    delete [] nnIdx;							// clean things up
    delete [] dists;
    annDeallocPt(queryPt);
  }
  


#ifdef _KEYPOINT_HPP_ 
//reprise de appariement descripteurs
  void NPPV(const std::vector<KeyPoint> & kptsdata, const std::vector<KeyPoint> & kptsquery, const unsigned int nbvois, std::vector<std::vector<unsigned int> > & idptsppv, std::vector<std::vector<double> > & distsppv, const int maxPtsAvisiter, const bool normalisedescripteur, const bool verbose){
    
    std::cout<<"Chargement des data..."<<std::flush;
    set_data(kptsdata);
    /*
    //juste pour test...
    std::vector<KeyPoint> data1,data2;
    for(int n=0;n<10;n++){data1.push_back(kptsdata[n]);}
    for(int n=10;n<kptsdata.size();n++){data2.push_back(kptsdata[n]);}
    set_size_dim(data1.size()+data2.size(),128);
    add_data(data1);
    add_data(data2);
    */

    std::cout<<"ok"<<std::endl;
    std::cout<<"Construction kdtree..."<<std::flush;
    construire_kdtree();
    std::cout<<"ok"<<std::endl;
    std::cout<<"Recherche ppv..."<<std::flush;
    nppv(kptsquery,nbvois,idptsppv,distsppv,maxPtsAvisiter,true);
    std::cout<<"ok"<<std::endl;
    
  }
#endif //_KEYPOINT_HPP_

protected:
  template<class TV>
  static const double Norme(const TV & pt, const unsigned char type_norme){
    double norme=0;
    if(type_norme==1){
      for(int n=0;n<pt.dim();n++){norme=norme+std::abs(pt[n]);}
    }else if(type_norme==2){
      for(int n=0;n<pt.dim();n++){norme=norme+std::pow(pt[n],2.0);}
      norme=std::sqrt(norme);
    }
    return norme;
  }

};



class RecherchePPV{
protected:
  unsigned int _NbPPV;
  unsigned int _NbmaxPtsAvisiter;
  bool _verbose;
public:
  //Constructeur / destructeur
  RecherchePPV():_NbPPV(1),_NbmaxPtsAvisiter(0),_verbose(true){}
  RecherchePPV(const unsigned int NbPPV, const unsigned int NbmaxPtsAvisiter, const bool verbose=true):_NbPPV(NbPPV),_NbmaxPtsAvisiter(NbmaxPtsAvisiter),_verbose(verbose){}
  virtual ~RecherchePPV(){}

  //Accesseurs 
  const unsigned int NbPPV() const{return _NbPPV;}
  const unsigned int NbmaxPtsAvisiter() const{return _NbmaxPtsAvisiter;}
  const bool verbose() const{return _verbose;}
  unsigned int & NbPPV() {return _NbPPV;}
  unsigned int & NbmaxPtsAvisiter() {return _NbmaxPtsAvisiter;}
  bool & verbose() {return _verbose;}

  //Fonctions
#ifdef _KEYPOINT_HPP_
  //recherche des nppv
  void nppv( KDTree & kdTree, const KeyPoint & kptquery, std::vector<unsigned int> & idptsppv, std::vector<double> & distsppv){
    int				nbnn				= _NbPPV;			// number of nearest neighbors
    double	       		eps				= 0;			// error bound
    
    if(kdTree.dim()!=(unsigned int)kptquery.dim()){
      std::cout<<"Dimension incompatible entre les points data et query"<<std::endl;
      return;
    }

    
    idptsppv.resize(_NbPPV);
    distsppv.resize(_NbPPV);
    

    //-------------------------------------------
    //Lecture des points query et recherche de leurs NN dans l'arbre
    //-------------------------------------------
    if(_verbose)std::cout<<"OK\nRequete..."<<std::flush;

    //Allocation
    ANNpoint		queryPt;				// query point
    queryPt = annAllocPt(kdTree.dim());	// allocate query point
    ANNidxArray		nnIdx;					// near neighbor indices
    ANNdistArray	        dists;					// near neighbor distances
    nnIdx = new ANNidx[nbnn];				// allocate near neigh indices
    dists = new ANNdist[nbnn];				// allocate near neighbor dists
    
    //if(normalisedescripteur)(*kptquery).normalise();
    for(unsigned int d=0;d<kdTree.dim();d++){
      queryPt[d]=kptquery[d];
    }
    
    annMaxPtsVisit(_NbmaxPtsAvisiter);//permet de limiter la recherche aux maxPtsAvisiter (par defaut 200) premiers points visités
    //kdTree->annkSearch(						// search
    kdTree.kdTree()->annkPriSearch(
			  queryPt,		// query point
			  nbnn,			// number of near neighbors
			  nnIdx,			// nearest neighbors (returned)
			  dists,			// distance (returned)
			  eps);			// error bound
    
    for(unsigned int k=0;k<_NbPPV;k++){
      idptsppv[k]=nnIdx[k];
      distsppv[k]=dists[k];
    }
    if(_verbose)std::cout<<"OK"<<std::endl;
    
    delete [] nnIdx;							// clean things up
    delete [] dists;
    annDeallocPt(queryPt);
  }
#endif //_KEYPOINT_HPP_

  


};



#endif
