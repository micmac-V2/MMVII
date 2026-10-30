#ifndef _FILTRAGEVOISINAGE_HPP_
#define _FILTRAGEVOISINAGE_HPP_

#include "ANN/ANN.h"					// ANN declarations

#include <list>
#include <vector>


//=============================================================
//Classe FiltreVoisinage
// Filtrage des appariements par consideration du voisinage
// en principe T=T2D<...>, mais peut aussi etre KeyPoint... désormais
//=============================================================
class FiltreVoisinage{
 protected:
  int _nbvoisins;//nombre de voisins a prendre en compte
  double _ratiovoisinsconserves;//proportion minimale de voisins conserves pour valider l'appariement

 public:
  //Constructeur / Destructeur
  FiltreVoisinage(const int nbvoisins=9,const double ratiovoisinsconserves=0.5):_nbvoisins(nbvoisins),_ratiovoisinsconserves(ratiovoisinsconserves){};
  virtual ~FiltreVoisinage(){}
  
  //Accesseurs
  const int nbvoisins() const{return _nbvoisins;}
  int & nbvoisins(){return _nbvoisins;}
  const double ratiovoisinsconserves() const{return _ratiovoisinsconserves;}
  double & ratiovoisinsconserves(){return _ratiovoisinsconserves;}

  //Fonctions
  template<class T>
    void filtre(std::vector<T> & pts1, std::vector<T> & pts2) const{
    FiltreVoisinage::FiltrageVoisinage(pts1,pts2,_ratiovoisinsconserves);
  }
  template<class T>
    void operator() (std::vector<T> & pts1, std::vector<T> & pts2) const{
    filtre(pts1,pts2);
  }
  template<class T>
    void filtre(std::vector<T> & pts1, std::vector<T> & pts2, const unsigned int nbiterations) const{
    for(unsigned int n=0;n<nbiterations;n++){
      unsigned int card=pts1.size();
      filtre(pts1,pts2);
      if(card==pts1.size()){break;}
    }
  }
  template<class T>
    void operator() (std::vector<T> & pts1, std::vector<T> & pts2, const unsigned int nbiterations) const{
    filtre(pts1,pts2,nbiterations);
  }

  template<class T>
  void getVoisinage(const std::vector<T> & pts, std::vector<std::vector<int> > & vois, const int nbvois=9, const bool verbose=true) const{
    FiltreVoisinage::Voisinage(pts,vois,nbvois,verbose);
  }


  //Fonctions statics

  //------------------------------------------------------------------
  //Calcule les nbvois plus proches voisins de chaque point du vecteur pts et les stocke dans vois
  //------------------------------------------------------------------
  template<class T>
    static void Voisinage(const std::vector<T> & pts, std::vector<std::vector<int> > & vois, const int nbvois=9, const bool verbose=true){
    if(pts.empty()){return;}
    //Parametres
    int dim=2; // dimension
    double eps=0; // error bound
    int maxPts=10000; // maximum number of data points
    //int maxPtsAvisiter = 200; // nombre maximal de points à visiter
    
    //on charge les points
    
    int			ndataPts;				// actual number of data points
    ANNpointArray		dataPts;				// data points
    ANNpoint		queryPt;				// query point
    ANNidxArray		nnIdx;					// near neighbor indices
    ANNdistArray	        dists;					// near neighbor distances
    ANNkd_tree*		kdTree;					// search structure
    
    queryPt = annAllocPt(dim);				// allocate query point
    //  dataPts = annAllocPts(maxPts, dim);			// allocate data points
    nnIdx=new ANNidx[nbvois];				// allocate near neigh indices
    dists=new ANNdist[nbvois];				// allocate near neighbor dists
    
    ndataPts=0;					        // read data points
    //Lecture des points data
    maxPts=pts.size();
    ndataPts=maxPts;
    dataPts=annAllocPts(maxPts, dim);			// allocate data points
    for(int i_data=0;i_data< ndataPts;i_data++) {					// read the data
      dataPts[i_data][0]=pts[i_data].x();
      dataPts[i_data][1]=pts[i_data].y();
    }
    
    //Construction du kd-tree
    if(verbose)std::cout<<"Construction de l'arbre..."<<std::flush;
    kdTree = new ANNkd_tree(				// build search structure
			    dataPts,			// the data points
			    ndataPts,			// number of points
			    dim);				// dimension of space
    
    //Lecture des points query et recherche de leurs NN dans l'arbre
    if(verbose)std::cout<<"Requetes..."<<std::flush;
    for(unsigned int i_query=0;i_query<pts.size();i_query++){
      queryPt[0]=pts[i_query].x();
      queryPt[1]=pts[i_query].y();
      
      //annMaxPtsVisit(maxPtsAvisiter);//permet de limiter la recherche aux maxPtsAvisiter (par defaut 200) premiers points visités
      //kdTree->annkSearch(						// search
      kdTree->annkPriSearch(
			    queryPt,		// query point
			    nbvois,			// number of near neighbors
			    nnIdx,			// nearest neighbors (returned)
			    dists,			// distance (returned)
			    eps);			// error bound
      
      //Recuperation de la requete
      std::list<int> ltmp;
      for(int i=0;i<nbvois;i++) {
	if(nnIdx[i]==(int)i_query){continue;}//chaque point est son propre voisin...
	ltmp.push_back(nnIdx[i]);//ltmp.push_back(nnIdx[i]);
      }
    //On peut reordonner (l'ordre actuel est celui des plus proches voisins...)
      ltmp.sort();
      std::vector<int> vtmp;for(std::list<int>::iterator itl=ltmp.begin();itl!=ltmp.end();++itl){vtmp.push_back((*itl));}
      vois.push_back(vtmp);
      
    }
    if(verbose)std::cout<<"OK"<<std::endl;
    
    delete [] nnIdx;							// clean things up
    delete [] dists;
    delete kdTree;
    annClose();									// done with ANN
  }
  


  //------------------------------------------------------------------
  //Effectue un filtrage des points faux en fonction de leurs voisins
  //------------------------------------------------------------------
  template<class T>
  static void FiltrageVoisinage(std::vector<T> & pts1, std::vector<T> & pts2, const double ratiovoisinsconserves=0.5){
    if(pts1.empty() || pts2.empty()){return;}
    std::vector<T> ptstmp1,ptstmp2;
    std::vector<std::vector<int> > vois1,vois2;
    Voisinage(pts1,vois1);
    Voisinage(pts2,vois2);
    for(unsigned int n=0;n<pts1.size();n++){
      //    SORTIEMESSAGE("\n"<<n<<" : "<<vois1[n][0]);for(unsigned int i=1;i<vois1[n].size();i++){SORTIEMESSAGE(" , "<<vois1[n][i]);}
      //    SORTIEMESSAGE("\n"<<n<<" : "<<vois2[n][0]);for(unsigned int i=1;i<vois2[n].size();i++){SORTIEMESSAGE(" , "<<vois2[n][i]);}
    }
    //Comparaison
    for(unsigned int n=0;n<pts1.size();n++){
      std::vector<int> v1,v2; v1=vois1[n]; v2=vois2[n];
      //Comparaison
      int nb=0;
      for(std::vector<int>::iterator i1=v1.begin();i1!=v1.end();++i1){
	for(std::vector<int>::iterator i2=v2.begin();i2!=v2.end();++i2){
	  if((*i1)==(*i2)){nb++;break;}
	}
      }
      if((double)(nb)>(ratiovoisinsconserves*(double)(v1.size()))){
	ptstmp1.push_back(pts1[n]);
	ptstmp2.push_back(pts2[n]);
      }
    }
    //  std::cout<<(pts1.size()-ptstmp1.size())<<" points supprimes sur "<<pts1.size()<<" points -> "<<((double)(pts1.size()-ptstmp1.size()))/((double)(pts1.size()))*100.0<<"%\n";
    pts1=ptstmp1;
    pts2=ptstmp2;
  }


};


#endif
