#ifndef _PARAMETRES_PRETRAITEMENT_DALLAGE_HPP_
#define _PARAMETRES_PRETRAITEMENT_DALLAGE_HPP_



//=========================================================
//ParametragePretraitement
//
//Classe de gestion des differents parametres lies a divers pretraitements de l'image : 
// - choix du canal ou de l'intensite
// - utilisation du gradient ? 
// - reetalement
//---------------------------------------------------------
class ParametragePretraitement{
protected :
  int    _canal;//si canal<0, intensite
  bool   _usegradient;//utilisation du gradient
  /*float*/ double  _maxim;
  /*float*/ double  _minim;//pas utilise pour le moment

public:
  ParametragePretraitement(){defaut();}
  virtual ~ParametragePretraitement(){};

  void defaut(){
    _canal          = 0 ;
    _usegradient    = false;
    _maxim          = -1;
    _minim          = 0;
  }

  //Parametres "radiometrie"
  const int canal() const{return _canal;}
  const bool usegradient() const{return _usegradient;}
  const /*float*/ double maxim() const{return _maxim;}
  const /*float*/ double minim() const{return _minim;}

  int & canal() {return _canal;}
  bool & usegradient() {return _usegradient;}
  /*float*/ double & maxim() {return _maxim;}
  /*float*/ double & minim() {return _minim;}

  
};


//=========================================================
//ParametrageDallage
//
//Classe de gestion des parametres parametres lies au decoupage en dalle de l'image a traiter
// - taille souhaitee (maximale) pour une dalle
// - taille minimale d'une dalle
// - marge optionnelle
//---------------------------------------------------------
class ParametrageDallage{
protected:
  //Parametrage du decoupage en dalles
  int _taille_souhaitee_dalle;//taille souhaitee (maximale) de la dalle
  int _taille_minimale_dalle;//taille minimale de la dalle
  int _marge;

public:
  ParametrageDallage(){defaut();}
  virtual ~ParametrageDallage(){};

  void defaut(){
    //Parametres decoupage en dalles
    _taille_souhaitee_dalle = 1024;
    _taille_minimale_dalle  = 200;
    _marge = 0;
  }

  //Parametres "decoupage"
  const int taille_souhaitee_dalle() const{return _taille_souhaitee_dalle;}
  const int taille_minimale_dalle() const{return _taille_minimale_dalle;}
  const int marge() const{return _marge;}

  int & taille_souhaitee_dalle(){return _taille_souhaitee_dalle;}
  int & taille_minimale_dalle(){return _taille_minimale_dalle;}
  int & marge() {return _marge;}
  
};


//=========================================================
//Dallage
//
// Choix des parametres pour decouper l'image en dalles de la maniere la plus reguliere possible
//----------------------------------------------
//On souhaite obtenir le decoupage le plus regulier possible de l'image en dalles :
//taille_image_a_decouper=taille_dalle*nb_dalles_identiques+taille_derniere_dalle
//avec taille_derniere_dalle>seuil
//---------------------------------------------------------
class Dallage{
protected:
  ParametrageDallage _param;
  bool _verbose;
public:
  Dallage(){}
  Dallage(const ParametrageDallage & param):_param(param){}
  //  Dallage(const ParametrageDallage & param, const double firstoctave):_param(param){
  //    _param.taille_souhaitee_dalle()=_param.taille_souhaitee_dalle()*(int)pow(2.0,firstoctave+1.0);
  //  }//potentiellement dangereux... -> memoire...
  virtual ~Dallage(){}

  const ParametrageDallage & param() const{return _param;}
  ParametrageDallage & param() {return _param;}
  const bool verbose() const{return _verbose;}
  bool & verbose(){return _verbose;}

  //==============================================
  // Choix des parametres pour decouper l'image en dalles de la maniere la plus reguliere possible
  //----------------------------------------------
  //On souhaite obtenir le decoupage le plus regulier possible de l'image en dalles :
  //taille_image_a_decouper=taille_dalle*nb_dalles_identiques+taille_derniere_dalle
  //avec taille_derniere_dalle>seuil
  
  static void ChoixMeilleursParametresDecoupe(const int taille_image_a_decouper, int & taille_dalle, int & nb_dalles_identiques, int & taille_derniere_dalle, const int seuil){
    if(taille_dalle<seuil){return; std::cout<<"Pas logique..."<<std::flush;}//pas logique...
    if(taille_image_a_decouper<taille_dalle){
      nb_dalles_identiques=1;
      taille_dalle=taille_image_a_decouper;
      taille_derniere_dalle=0; 
      return;
    }
    nb_dalles_identiques=(int)floor(((double)taille_image_a_decouper)/(double)taille_dalle);
    int dn=taille_image_a_decouper-taille_dalle*nb_dalles_identiques;
    if(dn==0){taille_derniere_dalle=0; return;}
    if(dn<seuil){
      taille_dalle=taille_dalle+(int)floor(((double)dn)/(double)nb_dalles_identiques);
      nb_dalles_identiques=(int)floor(((double)taille_image_a_decouper)/(double)taille_dalle);
      dn=taille_image_a_decouper-taille_dalle*nb_dalles_identiques;
      if(dn==0){taille_derniere_dalle=0; return;}
      nb_dalles_identiques=nb_dalles_identiques-1;
      taille_derniere_dalle=taille_dalle+dn;
    }else{
      taille_derniere_dalle=dn;
    }
  }

  static void ChoixMeilleursParametresDecoupe(const int taille_image_a_decouper, const int taille_dalle_souhaitee, const int seuil, std::vector<int> & tailles_dalles, const bool verbose=true){
    int taille_dalle=taille_dalle_souhaitee;
    int nb_dalles_identiques=0, taille_derniere_dalle=0;
    ChoixMeilleursParametresDecoupe(taille_image_a_decouper,taille_dalle,nb_dalles_identiques,taille_derniere_dalle,seuil);
    if(verbose)std::cout<< taille_dalle<<"*"<<nb_dalles_identiques<<"+"<<taille_derniere_dalle<<"="<<taille_dalle*nb_dalles_identiques+taille_derniere_dalle<<std::endl;  
    if(taille_derniere_dalle==0){
      tailles_dalles.resize(nb_dalles_identiques,taille_dalle);
    }else{
      tailles_dalles.resize(nb_dalles_identiques+1,taille_dalle);
      tailles_dalles[tailles_dalles.size()-1]=taille_derniere_dalle;
    } 
  }

  void ChoixMeilleursParametresDecoupe(const int taille_image_a_decouper, std::vector<int> & tailles_dalles) const{
    int seuil_taille_minimale=_param.taille_minimale_dalle();
    if(seuil_taille_minimale<=_param.marge()){
      seuil_taille_minimale=seuil_taille_minimale+_param.marge();
    }
    ChoixMeilleursParametresDecoupe(taille_image_a_decouper,_param.taille_souhaitee_dalle(),seuil_taille_minimale,tailles_dalles,_verbose);
  }
  
  //permet de tester la fonction ChoixMeilleursParametresDecoupe
  static void main_decoupage(int argc, char **argv){
    int N,n,nb,seuil,n_dernier;
    seuil=100;
    while(true){
      std::cin>>N>>n;
      ChoixMeilleursParametresDecoupe(N,n,nb,n_dernier,seuil);
      std::cout<<n<<"*"<<nb<<"+"<<n_dernier<<"="<<n*nb+n_dernier<<std::endl;
      std::vector<int> t;
      ChoixMeilleursParametresDecoupe(N,n,seuil,t);
      int somme=0;for(unsigned int i=0;i<t.size();i++){somme=somme+t[i];}
      std::cout<<somme<<std::endl;
    }
  }


};


#endif
