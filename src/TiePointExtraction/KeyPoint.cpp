#include "KeyPoint.h"

#include <iostream>



KeyPoint::KeyPoint(const int tailledesc):_desc(tailledesc,0){}

KeyPoint::KeyPoint(const double x, const double y, const int tailledesc):_desc(tailledesc,0){
  _x=x;
  _y=y;
}

KeyPoint::KeyPoint(const double x, const double y, const double sigma, const double orientation, const int tailledesc):_desc(tailledesc,0){
  _x=x;
  _y=y;
  _sigma=sigma;
  _orientation=orientation;
}

KeyPoint::KeyPoint(const double x, const double y, const double sigma, const double orientation, const bool signe, const int tailledesc):_desc(tailledesc,0){
  _x=x;
  _y=y;
  _sigma=sigma;
  _orientation=orientation;
  _signe=signe;
}

KeyPoint::KeyPoint(const double x, const double y, const double sigma, const double orientation, const std::vector<unsigned char> & desc):_desc(desc){
  _x=x;
  _y=y;
  _sigma=sigma;
  _orientation=orientation;
}

KeyPoint::KeyPoint(const double x, const double y, const double sigma, const double orientation, const bool signe, const std::vector<unsigned char> & desc):_desc(desc){
  _x=x;
  _y=y;
  _sigma=sigma;
  _orientation=orientation;
  _signe=signe;
}

KeyPoint::KeyPoint(const T2D<double> & xy, const int tailledesc){
  _x=xy.x();
  _y=xy.y();
  initdesc(tailledesc);
}

KeyPoint::KeyPoint(const T2D<double> & xy, const double sigma, const double orientation, const int tailledesc){
  _x=xy.x();
  _y=xy.y();
  _sigma=sigma;
  _orientation=orientation;
  initdesc(tailledesc);
}

KeyPoint::KeyPoint(const T2D<double> & xy, const double sigma, const double orientation, const bool signe, const int tailledesc){
  _x=xy.x();
  _y=xy.y();
  _sigma=sigma;
  _orientation=orientation;
  _signe=signe;
  initdesc(tailledesc);
}

KeyPoint::KeyPoint(const T2D<double> & xy, const double sigma, const double orientation, const std::vector<unsigned char> & desc):_desc(desc){
  _x=xy.x();
  _y=xy.y();
  _sigma=sigma;
  _orientation=orientation;
}

KeyPoint::KeyPoint(const T2D<double> & xy, const double sigma, const double orientation, const bool signe, const std::vector<unsigned char> & desc):_desc(desc){
  _x=xy.x();
  _y=xy.y();
  _sigma=sigma;
  _orientation=orientation;
  _signe=signe;
} 

const bool KeyPoint::operator == (const KeyPoint & kp){
  if(_x!=kp._x){return false;}
  if(_y!=kp._y){return false;}
  if(_sigma!=0 && kp._sigma!=0 &&_sigma!=kp._sigma){return false;}
  if(_orientation!=0 && kp._orientation!=0 &&_orientation!=kp._orientation){return false;}
  if(_desc!=kp._desc){return false;}//prevoir le cas ou on n'a qu'un point sans descripteur?
  return true;
}

const bool KeyPoint::operator != (const KeyPoint & kp){
  return !((*this)==kp);
}


inline const double KeyPoint::distanceeuclidiennecarree(const KeyPoint & kp1, const KeyPoint & kp2){
  double somme=0;
  int tailledesc=kp1.desc().size();
  if(tailledesc!=(int)kp2.desc().size()){return -1.0;}
  //  for(int i=0;i<tailledesc;i++){somme+=pow(kp1[i]-kp2[i],2.0);}
  KeyPoint::const_iterator i1,i2;
  i1=kp1.begin(); i2=kp2.begin();
  for(;i1!=kp1.end();++i1,++i2){somme+=pow((*i1)-(*i2),2.0);}
  return somme;
}
inline const double KeyPoint::distanceeuclidienne(const KeyPoint & kp1, const KeyPoint & kp2){
  return sqrt(distanceeuclidiennecarree(kp1,kp2));
}
inline const double KeyPoint::distance(const KeyPoint & kp1, const KeyPoint & kp2){
  return KeyPoint::distanceeuclidienne(kp1,kp2);
}
inline const double KeyPoint::distanceL1(const KeyPoint & kp1, const KeyPoint & kp2){
  double somme=0;
  int tailledesc=kp1.desc().size();
  if(tailledesc!=(int)kp2.desc().size()){return -1.0;}
  //  for(int i=0;i<tailledesc;i++){somme+=fabs(kp1[i]-kp2[i]);}
  KeyPoint::const_iterator i1,i2;
  i1=kp1.begin(); i2=kp2.begin();
  for(;i1!=kp1.end();++i1,++i2){somme+=fabs((double)(*i1)-(double)(*i2));}
  return sqrt(somme);
}

//permet de normaliser le vecteur descripteur
const std::vector<double> KeyPoint::normalise(const double newnorme) const {
  std::vector<double> descnormalise(_desc.size());
  std::vector<double>::iterator id=descnormalise.begin();
  KeyPoint::const_iterator it=_desc.begin();
  double n=newnorme/norme();
  //  for(it=_desc.begin();it!=_desc.end();++it){descnormalise.push_back(((double)(*it))/n);}
  for(;it!=_desc.end();++it,++id){(*id)=(((double)(*it))*n);}
  return descnormalise;
}

//on ne prend plus en compte la direction du gradient (approximation)
void KeyPoint::DemiDescripteur(const int nbins, const bool changetailledesc){//nbins=8 par defaut
  std::vector<double> desc;
  if(_desc.empty()){return;}
  if(changetailledesc){
    desc=std::vector<double>(dim()/2,0.0);
  }else{
    desc=std::vector<double>(dim(),0.0);
  }
  int moitie_nbins=nbins/2;
  double norme=0.0;
  double val=0.0;
  for(int n=0;n<dim()/nbins;n++){
    for(int k=0;k<moitie_nbins;k++){
      val=_desc[nbins*n+k]+_desc[nbins*n+k+moitie_nbins];
      if(changetailledesc){
      	desc[moitie_nbins*n+k]=val;
      }else{
	desc[nbins*n+k]=val;
      }
      norme=norme+pow(val,2);
    }
  }
  norme=255.0/sqrt(norme);

  if(changetailledesc){_desc.resize(dim()/2);}
  std::vector<double>::const_iterator id=desc.begin();
  iterator i=begin();
  for(;i!=end();++i,++id){(*i)=(unsigned char)(norme*(*id));}

}




//Fonctions de lecture et d'écriture des points d'intérêt et de leur descripteur dans des fichiers .key générés par l'exécutable de Lowe
void KeyPoint::lireKeyLowe(std::ifstream & f,const int tailledesc){
  f>>_y>>_x>>_sigma>>_orientation;
  if(_sigma<0){_sigma=-_sigma; _signe=false;}else{_signe=true;}
  for(int l=0;l<tailledesc;l++){int tmp;f>>tmp;_desc[l]=(unsigned char)tmp;}
}
void KeyPoint::lireKeyLowe(std::ifstream & f){lireKeyLowe(f,_desc.size());}

void KeyPoint::ecrireKeyLowe(std::ofstream & f) const{
  if(_signe){
    f<<_y<<" "<<_x<<" "<<_sigma<<" "<<_orientation<<std::endl;
  }else{
    f<<_y<<" "<<_x<<" "<<-_sigma<<" "<<_orientation<<std::endl;
  }
  for(int n=0;n<(int)_desc.size();n++){
    f<<" "<<(int)_desc[n];
    if(n%20==19){f<<std::endl;}
  }
  f<<"\n";
}

//Fonctions de lecture et d'écriture des points d'intérêt et de leur descripteur dans des fichiers .key générés par l'exécutable de l'INRIA
void KeyPoint::lireKeyINRIA(std::ifstream & f, const int tailledesc){
  double tmp;
  f>>_x>>_y>>_sigma>>_orientation>>tmp; //en fait les 3 derniers nombres sont les parametres permettant de definir une region affine...  a(x-_x)(x-_x)+2b(x-_x)(y-_y)+c(y-_y)(y-_y)=1
  if(_sigma<0){_sigma=-_sigma; _signe=false;}else{_signe=true;}
  for(int l=0;l<tailledesc;l++){int tmp;f>>tmp;_desc[l]=(unsigned char)tmp;}
}
void KeyPoint::lireKeyINRIA(std::ifstream & f){lireKeyLowe(f,_desc.size());}

void KeyPoint::ecrireKeyINRIA(std::ofstream & f) const{
  if(_signe){
    f<<_x<<" "<<_y<<" "<<_sigma<<" "<<_orientation<<" "<<"0"<<std::endl;
    //en fait les 3 derniers nombres doivent etre les parametres permettant de definir une region affine...  a(x-_x)(x-_x)+2b(x-_x)(y-_y)+c(y-_y)(y-_y)=1
  }else{
    f<<_y<<" "<<_x<<" "<<-_sigma<<" "<<_orientation<<" "<<"0"<<std::endl;
    //en fait les 3 derniers nombres doivent etre les parametres permettant de definir une region affine...  a(x-_x)(x-_x)+2b(x-_x)(y-_y)+c(y-_y)(y-_y)=1
  }
  for(int n=0;n<(int)dim();n++){
    f<<" "<<(int)_desc[n];
  }
  f<<std::endl;
}


void KeyPoint::ecrireKeyxml(std::ofstream & f) const{
  f<<"<KeypointN>"<<std::endl;
  f<<"\t<X>"<<_x<<"</X>"<<std::endl;
  f<<"\t<Y>"<<_y<<"</Y>"<<std::endl;
  if(_signe){
    f<<"\t<Scale>"<<_sigma<<"</Scale>"<<std::endl;
  }else{
    f<<"\t<Scale>"<<-_sigma<<"</Scale>"<<std::endl;
  }
  f<<"\t<Orientation>"<<_orientation<<"</Orientation>"<<std::endl;
  f<<"\t<Dim>"<<(int)dim()<<"</Dim>"<<std::endl;
  f<<"\t<Descriptor>"<<std::endl;
  for(int n=0;n<(int)dim();n++){
    f<<"\t\t<int>"<<(int)_desc[n]<<"</int>"<<std::endl;
  }
  f<<"\t</Descriptor>"<<std::endl;
  f<<"</KeypointN>"<<std::endl;
}


//Fonctions de lecture et d'écriture des points d'intérêt et de leur descripteur dans des fichiers binaires
void KeyPoint::lireBin(std::ifstream & f, const int tailledesc){
  float ftmp;
  f.read((char*)&ftmp,sizeof(float));_x=ftmp;
  f.read((char*)&ftmp,sizeof(float));_y=ftmp;
  f.read((char*)&ftmp,sizeof(float));_sigma=ftmp;
  f.read((char*)&ftmp,sizeof(float));_orientation=ftmp;
  if(_sigma<0){_sigma=-_sigma; _signe=false;}else{_signe=true;}
  //unsigned char uctmp;for(int l=0;l<tailledesc;l++){f.read((char*)&uctmp,sizeof(unsigned char));_desc[l]=uctmp;}
  //f.read((char*)&_desc[0],dim()*sizeof(unsigned char));
  for(iterator i=begin();i!=end();++i){f.read((char*)&(*i),sizeof(unsigned char));}
}
void KeyPoint::lireBin(std::ifstream & f){lireBin(f,_desc.size());}

void KeyPoint::ecrireBin(std::ofstream & f) const{
  float ftmp;
  ftmp=(float)_x;f.write((const char*)&ftmp,sizeof(float));
  ftmp=(float)_y;f.write((const char*)&ftmp,sizeof(float));
  ftmp=(float)_sigma; if(!_signe){ftmp=-ftmp;}  f.write((const char*)&ftmp,sizeof(float));
  ftmp=(float)_orientation;f.write((const char*)&ftmp,sizeof(float));
  //unsigned char uctmp;for(int n=0;n<(int)dim();n++){uctmp=(unsigned char)_desc[n];f.write((const char*)&uctmp,sizeof(unsigned char));}
  //f.write((const char*)&_desc[0],dim()*sizeof(unsigned char));
  for(const_iterator i=begin();i!=end();++i){f.write((const char*)&(*i),sizeof(unsigned char));}
}

//Fonctions de lecture et d'écriture des points d'intérêt et de leur descripteur dans des fichiers .key générés par l'exécutable siftpp
void KeyPoint::lireKeysiftpp(std::ifstream & f){
  f>>_x>>_y>>_orientation>>_sigma;
  if(_sigma<0){_sigma=-_sigma; _signe=false;}else{_signe=true;}
  for(unsigned int l=0;l<_desc.size();l++){int tmp;f>>tmp;_desc[l]=(unsigned char)tmp;}
}

void KeyPoint::ecrireKeysiftpp(std::ofstream & f) const{
  if(_signe){
    f<<_x<<" "<<_y<<" "<<_orientation<<" "<<_sigma;
  }else{
    f<<_x<<" "<<_y<<" "<<_orientation<<" "<<-_sigma;
  }
  for(unsigned int l=0;l<_desc.size();l++){f<<" "<<(int)_desc[l];}
  f<<"\n";
}

//Fonctions de lecture et d'écriture des points d'intérêt et de leur descripteur dans des fichiers .key générés par l'exécutable de Lowe
void KeyPoint::lireKeySURF(std::ifstream & f,const int tailledesc){
  double a,b,c;
  f>>_x>>_y>>a>>b>>c;
  double det = sqrt((a-c)*(a-c) + 4.0*b*b);
  double e1 = 0.5*(a+c + det);
  double e2 = 0.5*(a+c - det);
  double l1 = (1.0/sqrt(e1));
  double l2 = (1.0/sqrt(e2));
  double sc = sqrt( l1*l2 );
  _sigma=sc/2.5;
  double lap;
  f>>lap;
  if(lap>0){_signe=true;}else{_signe=false;}
  for(int l=0;l<tailledesc;l++){double tmp;f>>tmp;_desc[l]=(unsigned char)(tmp*255.0);}
}
void KeyPoint::lireKeySURF(std::ifstream & f){lireKeyLowe(f,_desc.size());}

/*
void KeyPoint::ecrireKeySURF(std::ofstream & f) const{
  f<<_y<<" "<<_x<<" "<<_sigma<<" "<<_orientation<<"\n";
  for(int n=0;n<(int)_desc.size();n++){
    f<<" "<<(int)_desc[n];
    if(n%20==19){f<<"\n";}
  }
  f<<"\n";
  }*/




//Fonctions de lecture et d'écriture du point d'intérêt et de son descripteur dans des fichiers .loc et .pts
void KeyPoint::ecrireLoc(std::ofstream & f) const{f<<_x<<"  "<<_y<</*"  "<<*/"\n";}//inversion x y corrigee

void KeyPoint::ecrireDesc(std::ofstream & f) const{
  //  double N=norme();
  for(std::vector<unsigned char>::const_iterator it=_desc.begin();it!=_desc.end();++it){
    f<<(int)(*it)<<"  ";
    //    f<<(((double)(*it))/N)<<"  ";//modif passage unsigned char
  }
  f<<"\n";
}

void KeyPoint::ecrireLocDesc(std::ofstream & floc, std::ofstream & fdesc) const{ecrireLoc(floc);ecrireDesc(fdesc);}
void KeyPoint::lireLocDesc(std::ifstream & floc, std::ifstream & fdesc, const int tailledesc){
  floc>>_x>>_y;//x et y inverses corrige!
  _desc=std::vector<unsigned char>(tailledesc);
  int tmp;
  for(int n=0;n<tailledesc;n++){
    fdesc>>tmp;_desc[n]=(unsigned char)tmp;
  }
}
void KeyPoint::lireLocDesc(std::ifstream & floc, std::ifstream & fdesc){lireLocDesc(floc,fdesc,tailledesc());}

const double KeyPoint::norme() const{   
  double somme=0;
  std::vector<unsigned char>::const_iterator it;
  for(it=_desc.begin();it!=_desc.end();++it){somme+=pow((double)(*it),2.0);}
  return sqrt(somme);
}//calcul de la norme du vecteur descipteur SIFT


//-------------------------------------------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Namespace ManipKeyPoint de manipulation (lecture-écriture, ...) de points d'intérêts (objets KeyPoint)
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------------------------------------


const int OutilsIOKeyPoint::card_KeyLowe(const stringpp & adresse_key){
  std::ifstream fk(adresse_key.c_str(),std::ios::in);
  if(!fk.good()){return 0;}
  int nbpts;
  fk>>nbpts;
  fk.close();
  return nbpts;
}

const int OutilsIOKeyPoint::tailledesc_KeyLowe(const stringpp & adresse_key){
  std::ifstream fk(adresse_key.c_str(),std::ios::in);
  if(!fk.good()){return 0;}
  int tailledesc;
  fk>>tailledesc>>tailledesc;
  fk.close();
  return tailledesc;
}

void OutilsIOKeyPoint::entete_KeyLowe(const stringpp & adresse_key, int & nbpts, int & tailledesc){
  std::ifstream fk(adresse_key.c_str(),std::ios::in);
  if(!fk.good()){return;}
  fk>>nbpts>>tailledesc;
  fk.close();
}

const int OutilsIOKeyPoint::card_KeyBin(const stringpp & adresse_key){
  std::ifstream fk(adresse_key.c_str(),std::ios::in | std::ios::binary);
  return card_KeyBin(fk);
}
const int OutilsIOKeyPoint::card_KeyBin(std::ifstream & f_key){
  if(!f_key.good()){return 0;}
  int nbpts;
  f_key.read((char*)&nbpts,sizeof(unsigned int));
  return nbpts;
}
const int OutilsIOKeyPoint::tailledesc_KeyBin(const stringpp & adresse_key){
  std::ifstream fk(adresse_key.c_str(),std::ios::in | std::ios::binary);
  return tailledesc_KeyBin(fk);
}
const int OutilsIOKeyPoint::tailledesc_KeyBin(std::ifstream & f_key){
  if(!f_key.good()){return 0;}
  int tailledesc;
  f_key.read((char*)&tailledesc,sizeof(unsigned int));
  f_key.read((char*)&tailledesc,sizeof(unsigned int));
  return tailledesc;
}
void OutilsIOKeyPoint::entete_KeyBin(const stringpp & adresse_key, int & nbpts, int & tailledesc){
  std::ifstream fk(adresse_key.c_str(),std::ios::in | std::ios::binary);
  entete_KeyBin(fk,nbpts,tailledesc);
}
void OutilsIOKeyPoint::entete_KeyBin(std::ifstream & f_key, int & nbpts, int & tailledesc){
  if(!f_key.good()){return;}
  f_key.read((char*)&nbpts,sizeof(unsigned int));
  f_key.read((char*)&tailledesc,sizeof(unsigned int));
}


const int OutilsIOKeyPoint::card_KeyINRIA(const stringpp & adresse_key){
  std::ifstream fk(adresse_key.c_str(),std::ios::in);
  if(!fk.good()){return 0;}
  int nbpts;
  fk>>nbpts>>nbpts;
  fk.close();
  return nbpts;
}

const int OutilsIOKeyPoint::tailledesc_KeyINRIA(const stringpp & adresse_key){
  std::ifstream fk(adresse_key.c_str(),std::ios::in);
  if(!fk.good()){return 0;}
  int tailledesc;
  fk>>tailledesc;
  fk.close();
  return tailledesc;
}

void OutilsIOKeyPoint::entete_KeyINRIA(const stringpp & adresse_key, int & nbpts, int & tailledesc){
  std::ifstream fk(adresse_key.c_str(),std::ios::in);
  if(!fk.good()){return;}
  fk>>tailledesc>>nbpts;
  fk.close();
}



const int OutilsIOKeyPoint::card_Loc(const stringpp & adresse_loc){
  std::ifstream floc(adresse_loc.c_str(),std::ios::in);
  if(!floc.good()){return 0;}
  int nbpts=0;
  double tmp;
  while(true){
    floc>>tmp>>tmp;
    if(floc.eof()){break;}
    nbpts++;
  }
  return nbpts;
}

const int OutilsIOKeyPoint::card_Desc(const stringpp & adresse_desc, const int tailledesc){
  std::ifstream fdesc(adresse_desc.c_str(),std::ios::in);
  if(!fdesc.good()){return 0;}
  int nbpts=0;
  /*double*/int tmp;
  while(true){
    for(int k=0;k<tailledesc;k++){fdesc>>tmp;}
    if(fdesc.eof()){break;}
    nbpts++;
  }
  return nbpts;
}
 
const int OutilsIOKeyPoint::card_LocDesc(const stringpp & adresse_loc, const stringpp & adresse_desc, const int tailledesc){
  std::ifstream floc(adresse_loc.c_str(),std::ios::in);
  std::ifstream fdesc(adresse_desc.c_str(),std::ios::in);
  if(!(floc.good() && fdesc.good())){return 0;}
  int nbpts=0;
  double tmp;
  while(true){
    floc>>tmp>>tmp;
    for(int k=0;k<tailledesc;k++){fdesc>>tmp;}
    if(floc.eof() && fdesc.eof()){break;}
    else if(floc.eof() || fdesc.eof()){stringpp err("Les fichiers '.loc' et '.pts' n'ont pas la meme taille..."); std::cout<<err; throw err; return -1;}
    nbpts++;
  }
  return nbpts;
} 

const int OutilsIOKeyPoint::card_LocDesc(const stringpp & adresse, const int tailledesc){
  return OutilsIOKeyPoint::card_LocDesc(adresse.newSetExtension(".loc"),adresse.newSetExtension(".desc"),tailledesc);
}


//---------------------------------------------------
// Traduction d'un fichier .key (généré par l'exécutable de Lowe) en un fichier .loc et en un fichier .pts
//---------------------------------------------------
void OutilsIOKeyPoint::KeyLowe2LocDesc(const stringpp & adresse_key, const stringpp & adresse_loc, const stringpp & adresse_desc){
  int nbpts=0;
  OutilsIOKeyPoint::KeyLowe2LocDesc_card(adresse_key, adresse_loc, adresse_desc, nbpts);
}

void OutilsIOKeyPoint::KeyLowe2LocDesc(const stringpp & adresse){
  OutilsIOKeyPoint::KeyLowe2LocDesc(adresse.newSetExtension(".key"),adresse.newSetExtension(".desc"),adresse.newSetExtension(".pts"));
}

void OutilsIOKeyPoint::KeyLowe2LocDesc_card(const stringpp & adresse_key, const stringpp & adresse_loc, const stringpp & adresse_desc, int & nbpts){
  int tailledesc=0;
  OutilsIOKeyPoint::KeyLowe2LocDesc_card(adresse_key, adresse_loc, adresse_desc, nbpts, tailledesc);
}

void OutilsIOKeyPoint::KeyLowe2LocDesc_card(const stringpp & adresse, int & nbpts){
  OutilsIOKeyPoint::KeyLowe2LocDesc_card(adresse.newSetExtension(".key"),adresse.newSetExtension(".loc"),adresse.newSetExtension(".desc"),nbpts);
}


void OutilsIOKeyPoint::KeyLowe2LocDesc_card(const stringpp & adresse_key, const stringpp & adresse_loc, const stringpp & adresse_desc, int & nbpts, int & tailledesc){
  std::ifstream fkey(adresse_key.c_str(),std::ios::in);
  std::ofstream floc(adresse_loc.c_str(),std::ios::out);
  std::ofstream fdesc(adresse_desc.c_str(),std::ios::out);

  if(!fkey.good()){return;}
  nbpts=0;
  tailledesc=0;
  fkey>>nbpts;
  fkey>>tailledesc;
  for(int n=0;n<nbpts;n++){
    if(fkey.eof()){std::cout<<"Pb:Il y a moins de points que prevu.\n"; return;}
    //Mieux vaut utiliser la classe KeyPoint (cf normalisation du vecteur desc dans le fichier .pts ...)
    KeyPoint kptmp(tailledesc);
    kptmp.lireKeyLowe(fkey);
    kptmp.ecrireLocDesc(floc,fdesc);
  }
  floc.close();fdesc.close();fkey.close();
}

void OutilsIOKeyPoint::KeyLowe2LocDesc_card(const stringpp & adresse, int & nbpts, int & tailledesc){
  OutilsIOKeyPoint::KeyLowe2LocDesc_card(adresse.newSetExtension(".key"),adresse.newSetExtension(".loc"),adresse.newSetExtension(".desc"),nbpts,tailledesc);
}


//--------------------------------------------------------
//Traduction d'un fichier .loc et d'un fichier .pts en un fichier .key (du type de ceux générés par l'exécutable de Lowe)
//--------------------------------------------------------
void OutilsIOKeyPoint::LocDesc2KeyLowe(const stringpp & adresse_loc, const stringpp & adresse_desc, const stringpp & adresse_key, const int tailledesc){
  int nbpts=0;
  LocDesc2KeyLowe_card(adresse_loc, adresse_desc, adresse_key, nbpts, tailledesc);
}

void OutilsIOKeyPoint::LocDesc2KeyLowe(const stringpp & adresse, const int tailledesc){
  OutilsIOKeyPoint::LocDesc2KeyLowe(adresse.newSetExtension(".loc"), adresse.newSetExtension(".desc"), adresse.newSetExtension(".key"), tailledesc);
}

void OutilsIOKeyPoint::LocDesc2KeyLowe_card(const stringpp & adresse_loc, const stringpp & adresse_desc, const stringpp & adresse_key, int & nbpts, const int tailledesc){
  std::ifstream floc1(adresse_loc.c_str(), std::ios::in);
  std::ifstream fdesc1(adresse_desc.c_str(), std::ios::in);

  if(!(floc1.good() && fdesc1.good())){return;}

  nbpts=0;
  
  bool doubleverification=false;
  double tmp;
  while(true){
    floc1>>tmp>>tmp;
    if(doubleverification){for(int n=0;n<tailledesc;n++){fdesc1>>tmp;}}
    if(doubleverification){
      if(floc1.eof() && fdesc1.eof()){break;}
      else if(floc1.eof() || fdesc1.eof()){stringpp err("Les fichiers '.loc' et '.desc' n'ont pas la meme taille...");std::cout<<err;throw err;return;}
    }else{
      if(floc1.eof()){break;}
    }
    nbpts++;
  }
  
  //  fdesc.seekg(0); //pour revenir au début de fdesc
  //  floc.seekg((std::streamoff)(0));

  std::ofstream fkey(adresse_key.c_str(), std::ios::out);
  std::ifstream floc(adresse_loc.c_str(), std::ios::in);
  std::ifstream fdesc(adresse_desc.c_str(), std::ios::in);

  if(!(floc.good() && fdesc.good())){return;}

  fkey<<nbpts<<" "<<tailledesc<<"\n";
  KeyPoint kptmp(tailledesc);
  for(int n=0;n<nbpts;n++){
    kptmp.lireLocDesc(floc,fdesc);
    kptmp.ecrireKeyLowe(fkey);
  }
  
  fkey.close();floc.close();fdesc.close();
  //*/

  /*//en une seule fois
  std::ifstream floc(adresse_loc, std::ios::in);
  std::ifstream fdesc(adresse_desc, std::ios::in);
  std::ofstream fkey(adresse_key, std::ios::out);
  fkey<<"     ";
  fkey<<"     ";//fkey<<"\n\n\n\n\n";
  fkey<<" "<<tailledesc<<"\n";
  int nbpts=0;
  while(true){
    KeyPoint kptmp(tailledesc);
    kptmp.lireLocPts(floc,fdesc);
    if(fdesc.eof() && floc.eof()){break;}
    else if(fdesc.eof() || floc.eof()){return;}//erreur les fichiers n'ont pas la meme taille...
    kptmp.ecrireKeyLowe(fkey);
    nbpts++;
  } 
  //  floc.seekg(0,std::ios::beg); //pour revenir au début de floc
  //  fdesc.seekg(0,std::ios::beg); //pour revenir au début de fdesc
  //  fkey.seekp(0,std::ios::beg);
  fkey.seekp(0);
  fkey<<nbpts;
  fdesc.close();floc.close();fkey.close();
  //*/
}

void OutilsIOKeyPoint::LocDesc2KeyLowe_card(const stringpp & adresse, int & nbpts, const int tailledesc){
  OutilsIOKeyPoint::LocDesc2KeyLowe_card(adresse.newSetExtension(".loc"), adresse.newSetExtension(".desc"), adresse.newSetExtension(".key"), nbpts, tailledesc);
}


//-----------------------------------------------------------------
//Chargement d'un fichier .key (généré par l'exécutable de Lowe). On obtient en sortie un vector de points d'intérêt (objets KeyPoint)
//-----------------------------------------------------------------
const bool OutilsIOKeyPoint::ChargerKeyLowe(const stringpp & adresse_key, std::vector<KeyPoint> & kp){
  std::ifstream Fkey(adresse_key.c_str(),std::ios::in);
  if(!Fkey.good()){std::cout<<"Probleme en lecture du fichier "<<adresse_key<<std::endl; return false;}
  int nbpts,tailledesc;
  Fkey>>nbpts;
  kp.resize(nbpts);
  Fkey>>tailledesc;
  for(int n=0;n<nbpts;n++){
    if(Fkey.eof()){std::cout<<"Pb:Il y a moins de points que prevu.\n";return false;}
    KeyPoint kptmp(tailledesc);
    kptmp.lireKeyLowe(Fkey);
    kp[n]=kptmp;
  }
  Fkey.close();
  return true;
}

const bool OutilsIOKeyPoint::ChargerKeyLowe(const stringpp & adresse_key, std::vector<T2D<double> > & kp){
  std::ifstream Fkey(adresse_key.c_str(),std::ios::in);
  if(!Fkey.good()){std::cout<<"Probleme en lecture du fichier "<<adresse_key<<std::endl; return false;}
  int nbpts,tailledesc;
  double tmp;
  Fkey>>nbpts;
  kp.resize(nbpts);
  Fkey>>tailledesc;
  for(int n=0;n<nbpts;n++){
    if(Fkey.eof()){std::cout<<"Pb:Il y a moins de points que prevu.\n"; return false;}
    double x,y;
    Fkey>>y>>x;//ATTENTION a l'inversion!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    kp[n]=T2D<double> (x,y);//kp.push_back(T2D<double> (x,y));
    Fkey>>tmp>>tmp;
    for(int k=0;k<tailledesc;k++){
      if(Fkey.eof()){std::cout<<"Pb:Il y a moins de points que prevu.\n"; return false;}
      Fkey>>tmp;
    }
  }
  return true;
}

const bool OutilsIOKeyPoint::ChargerKeyINRIA(const stringpp & adresse_key, std::vector<KeyPoint> & kp){
  std::ifstream Fkey(adresse_key.c_str(),std::ios::in);
  if(!Fkey.good()){std::cout<<"Probleme en lecture du fichier "<<adresse_key<<std::endl; return false;}
  int nbpts,tailledesc;
  Fkey>>tailledesc;
  Fkey>>nbpts;
  kp.resize(nbpts);
  for(int n=0;n<nbpts;n++){
    if(Fkey.eof()){std::cout<<"Pb:Il y a moins de points que prevu.\n"; return false;}
    KeyPoint kptmp(tailledesc);
    kptmp.lireKeyINRIA(Fkey);
    kp[n]=kptmp;//kp.push_back(kptmp);
  }
  Fkey.close();
  return true;
}

const bool OutilsIOKeyPoint::ChargerKeyINRIA(const stringpp & adresse_key, std::vector<T2D<double> > & kp){
  std::ifstream Fkey(adresse_key.c_str(),std::ios::in);
  if(!Fkey.good()){std::cout<<"Probleme en lecture du fichier "<<adresse_key<<std::endl; return false;}
  int nbpts,tailledesc;
  double tmp;
  Fkey>>tailledesc;
  Fkey>>nbpts;
  kp.resize(nbpts);
  for(int n=0;n<nbpts;n++){
    if(Fkey.eof()){std::cout<<"Pb:Il y a moins de points que prevu.\n" ;return false;}
    double x,y;
    Fkey>>x>>y;//ATTENTION a l'inversion!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    kp[n]=T2D<double>(x,y);//kp.push_back(T2D<double> (x,y));
    Fkey>>tmp>>tmp>>tmp;
    for(int k=0;k<tailledesc;k++){
      if(Fkey.eof()){std::cout<<"Pb:Il y a moins de points que prevu.\n"; return false;}
      Fkey>>tmp;
    }
  }
  return true;
}


const bool OutilsIOKeyPoint::ChargerKeyBin(const stringpp & adresse_key, std::vector<KeyPoint> & kp){
  std::ifstream Fkey(adresse_key.c_str(),std::ios::in | std::ios::binary);
  if(!Fkey.good()){std::cout<<"Probleme en lecture du fichier "<<adresse_key<<std::endl; return false;}
  unsigned int nbpts, tailledesc;
  Fkey.read((char*)&nbpts,sizeof(unsigned int));
  Fkey.read((char*)&tailledesc,sizeof(unsigned int));
  kp.resize(nbpts);
  for(unsigned int n=0;n<nbpts;n++){
    if(Fkey.eof()){std::cout<<"Pb:Il y a moins de points que prevu.\n"; return false;}
    KeyPoint kptmp(tailledesc);
    kptmp.lireBin(Fkey);
    kp[n]=kptmp;//kp.push_back(kptmp);
  }
  Fkey.close();
  return true;
}

const bool OutilsIOKeyPoint::ChargerKeyBin(const stringpp & adresse_key, std::vector<KeyPoint> & kp, const double sigmamin, const double sigmamax){
  std::ifstream Fkey(adresse_key.c_str(),std::ios::in | std::ios::binary);
  if(!Fkey.good()){std::cout<<"Probleme en lecture du fichier "<<adresse_key<<std::endl; return false;}
  unsigned int nbpts, tailledesc;
  Fkey.read((char*)&nbpts,sizeof(unsigned int));
  Fkey.read((char*)&tailledesc,sizeof(unsigned int));
  kp.resize(nbpts);
  unsigned int nbptsgardes=0;
  for(unsigned int n=0;n<nbpts;n++){
    if(Fkey.eof()){std::cout<<"Pb:Il y a moins de points que prevu.\n"; return false;}
    KeyPoint kptmp(tailledesc);
    kptmp.lireBin(Fkey);
    if(fabs(kptmp.sigma())<sigmamin){continue;}
    if(fabs(kptmp.sigma())>sigmamax){continue;}
    kp[nbptsgardes]=kptmp;//kp.push_back(kptmp);
    nbptsgardes++;
  }
  kp.resize(nbptsgardes);
  Fkey.close();
  return true;
}


const bool OutilsIOKeyPoint::ChargerKeyBin(const stringpp & adresse_key, std::vector<T2D<double> > & kp){
  std::ifstream Fkey(adresse_key.c_str(),std::ios::in | std::ios::binary);
  if(!Fkey.good()){std::cout<<"Probleme en lecture du fichier "<<adresse_key<<std::endl; return false;}
  unsigned int nbpts, tailledesc;
  Fkey.read((char*)&nbpts,sizeof(unsigned int));
  Fkey.read((char*)&tailledesc,sizeof(unsigned int));
  kp.resize(nbpts);
  for(unsigned int n=0;n<nbpts;n++){
    if(Fkey.eof()){std::cout<<"Pb:Il y a moins de points que prevu.\n"; return false;}
    KeyPoint kptmp(tailledesc);
    kptmp.lireBin(Fkey);
    kp[n]=T2D<double> (kptmp.x(),kptmp.y());//kp.push_back(T2D<double> (kptmp.x(),kptmp.y()));
  }
  return true;
}



const bool OutilsIOKeyPoint::ChargerLoc(const stringpp & adresse_loc, std::vector<T2D<double> > & kp){
  std::ifstream Floc(adresse_loc.c_str(),std::ios::in);
  if(!Floc.good()){std::cout<<"Probleme en lecture du fichier "<<adresse_loc<<std::endl; return false;}
  while(true){
    double x,y;
    Floc>>x>>y;
    if(Floc.eof()){break;}
    kp.push_back(T2D<double>(x,y));
  }
  return true;
}


//------------------------------------------------------------------------
//Chargement d'un fichier .loc et d'un fichier .pts . On obtient en sortie un vector de points d'intérêt (objets KeyPoint)
//------------------------------------------------------------------------
const bool OutilsIOKeyPoint::ChargerLocDesc(const stringpp & adresse_loc, const stringpp & adresse_desc, std::vector<KeyPoint> & kp, const int tailledesc){
  std::ifstream floc(adresse_loc.c_str(),std::ios::in);
  std::ifstream fdesc(adresse_desc.c_str(),std::ios::in);
  if(!(floc.good() && fdesc.good())){std::cout<<"Probleme en lecture des fichiers "<<adresse_loc<<" et "<<adresse_desc<<std::endl; return false;}
  while(true){
    KeyPoint kptmp(tailledesc);
    kptmp.lireLocDesc(floc,fdesc);
    if(fdesc.eof() && floc.eof()){break;}
    if(fdesc.eof() || floc.eof()){stringpp err("Les fichiers '.loc' et '.desc' n'ont pas la meme taille..."); std::cout<<err; throw err;}
    kp.push_back(kptmp);
  }
  floc.close();fdesc.close();
  return true;
}

const bool OutilsIOKeyPoint::ChargerLocDesc(const stringpp & adresse, std::vector<KeyPoint> & kp, const int tailledesc){
  return OutilsIOKeyPoint::ChargerLocDesc(adresse.newSetExtension(".loc"), adresse.newSetExtension(".desc"), kp, tailledesc);
}


//-------------------------------------------------------------------------
//Ecriture d'une liste de points d'intérêt et leur descripteur (objets KeyPoint) sous forme d'un fichier .key (du type de ceux générés par l'exécutable de Lowe)
//-------------------------------------------------------------------------
const bool OutilsIOKeyPoint::EcrireKeyLowe(const stringpp & adresse_key, const std::vector<KeyPoint> & kp){
  std::ofstream f(adresse_key.c_str(), std::ios::out);
  if(!f.good()){std::cout<<"Probleme en ecriture du fichier "<<adresse_key<<std::endl; return false;}
  f<<kp.size()<<" "<<kp[0].tailledesc()<<"\n";
  std::vector<KeyPoint>::const_iterator it;
  for(it=kp.begin();it!=kp.end();++it){(*it).ecrireKeyLowe(f);}
  f.close();
  return true;
}

const bool OutilsIOKeyPoint::EcrireKeyINRIA(const stringpp & adresse_key, const std::vector<KeyPoint> & kp){
  std::ofstream f(adresse_key.c_str(), std::ios::out);
  if(!f.good()){std::cout<<"Probleme en ecriture du fichier "<<adresse_key<<std::endl; return false;}
  f<<kp[0].tailledesc()<<"\n"<<kp.size()<<"\n";
  std::vector<KeyPoint>::const_iterator it;
  for(it=kp.begin();it!=kp.end();++it){(*it).ecrireKeyINRIA(f);}
  f.close();
  return true;
}

const bool OutilsIOKeyPoint::EcrireKeyBin(const stringpp & adresse_key, const std::vector<KeyPoint> & kp){
  std::ofstream f(adresse_key.c_str(), std::ios::out | std::ios::binary);
  if(!f.good()){std::cout<<"Probleme en ecriture du fichier "<<adresse_key<<std::endl; return false;}
  unsigned int nbpts=kp.size();
  unsigned int tdesc=kp[0].tailledesc();
  f.write((const char*)&nbpts,sizeof(unsigned int));
  f.write((const char*)&tdesc,sizeof(unsigned int));
  std::vector<KeyPoint>::const_iterator it;
  for(it=kp.begin();it!=kp.end();++it){(*it).ecrireBin(f);}
  f.close();
  return true;
}

const bool OutilsIOKeyPoint::EcrireKeyxml(const stringpp & adresse_key, const std::vector<KeyPoint> & kp){
  std::ofstream f(adresse_key.c_str(), std::ios::out);
  if(!f.good()){std::cout<<"Probleme en ecriture du fichier "<<adresse_key<<std::endl; return false;}
  std::vector<KeyPoint>::const_iterator it;
  for(it=kp.begin();it!=kp.end();++it){(*it).ecrireKeyxml(f);}
  f.close();
  return true;
}

//-------------------------------------------------------------------------
//Ecriture d'une liste de points d'intérêt (uniquement leurs coordonnées...) donnée sous la forme d'un vector d'objets KeyPoint dans un fichier .loc
//-------------------------------------------------------------------------
const bool OutilsIOKeyPoint::EcrireLoc(const stringpp & adresse_loc, const std::vector<KeyPoint> & kp){
  std::ofstream loc(adresse_loc.c_str(), std::ios::out);
  if(!loc.good()){std::cout<<"Probleme en ecriture du fichier "<<adresse_loc<<std::endl; return false;}
  std::vector<KeyPoint>::const_iterator it;
  for(it=kp.begin();it!=kp.end();++it){(*it).ecrireLoc(loc);}
  loc.close();
  return true;
}

//-------------------------------------------------------------------------
//Ecriture d'une liste de points d'intérêt (uniquement leur descripteur...) donnée sous la forme d'un vector d'objets KeyPoint dans un fichier .pts
//-------------------------------------------------------------------------
const bool OutilsIOKeyPoint::EcrireDesc(const stringpp & adresse_desc, const std::vector<KeyPoint> & kp){
  std::ofstream fdesc(adresse_desc.c_str(), std::ios::out);
  if(!fdesc.good()){std::cout<<"Probleme en ecriture du fichier "<<adresse_desc<<std::endl; return false;}
  std::vector<KeyPoint>::const_iterator it;
  for(it=kp.begin();it!=kp.end();++it){(*it).ecrireDesc(fdesc);}
  fdesc.close();
  return true;
}

//---------------------------------------------------------------------------
//Ecriture d'une liste de points d'intérêt donnés sous la forme d'un vector d'objets KeyPoint dans un fichier .loc (pour leurs coordonnées) et dans un fichier .pts (pour leur descripteur)
//---------------------------------------------------------------------------
const bool OutilsIOKeyPoint::EcrireLocDesc(const stringpp & adresse_loc, const stringpp & adresse_desc, const std::vector<KeyPoint> & kp){
  std::ofstream floc(adresse_loc.c_str(), std::ios::out);
  std::ofstream fdesc(adresse_desc.c_str(), std::ios::out);
  if(!floc.good() || !fdesc.good()){std::cout<<"Probleme en ecriture des fichiers "<<adresse_loc<<" et "<<adresse_desc<<std::endl; return false;}
  std::vector<KeyPoint>::const_iterator it;
  for(it=kp.begin();it!=kp.end();++it){(*it).ecrireDesc(fdesc);(*it).ecrireLoc(floc);}
  floc.close();fdesc.close();
  return true;
}

const bool OutilsIOKeyPoint::EcrireLocDesc(const stringpp & adresse, const std::vector<KeyPoint> & kp){
  return OutilsIOKeyPoint::EcrireLocDesc(adresse.newSetExtension(".loc"), adresse.newSetExtension(".desc"), kp);
}


//---------------------------------------------------------------------------
//Calcul de la distance entre les descripteurs de deux points d'intérêt
//---------------------------------------------------------------------------
const double OutilsKeyPoint::distanceeuclidienne(const KeyPoint & kp1, const KeyPoint & kp2){
  return KeyPoint::distanceeuclidienne(kp1,kp2);
}
const double OutilsKeyPoint::distance(const KeyPoint & kp1, const KeyPoint & kp2){
  return KeyPoint::distance(kp1,kp2);
}
const double OutilsKeyPoint::distanceL1(const KeyPoint & kp1, const KeyPoint & kp2){
  return KeyPoint::distanceL1(kp1,kp2);
}


//--------------------------------------------------------------
//Traduction  d'un fichier siftpp en un fichier sift de Lowe
//--------------------------------------------------------------
void OutilsIOKeyPoint::Siftpp2KeyLowe(const stringpp & adresse_siftpp, const stringpp & adresse_keylowe){
  std::vector<KeyPoint> kp;
  std::ifstream ksiftpp(adresse_siftpp.c_str(),std::ios::in);
  if(!ksiftpp.good()){return;}
  std::ofstream kplowe(adresse_keylowe.c_str(),std::ios::out);
  while(true){
    KeyPoint kptmp;
    kptmp.lireKeysiftpp(ksiftpp);
    if(ksiftpp.eof()){break;}
    kp.push_back(kptmp);
  }
  kplowe<<kp.size()<<" 128\n";
  for(unsigned int i=0;i<kp.size();i++){
    kp[i].ecrireKeyLowe(kplowe);
  }

  ksiftpp.close();
  kplowe.close();
}


//--------------------------------------------------------------
//Traduction  d'un fichier sift de Lowe <->fichier binaire
//--------------------------------------------------------------  
void OutilsIOKeyPoint::KeyLowe2KeyBin(const stringpp & adresse_keylowe, const stringpp & adresse_keybin){
  std::ifstream kl(adresse_keylowe.c_str(),std::ios::in);
  std::ofstream kb(adresse_keybin.c_str(),std::ios::out | std::ios::binary);
  if(!kl.good()){return;}
  unsigned int nbpts, tdesc;
  kl>>nbpts>>tdesc;
  kb.write((const char*)&nbpts,sizeof(unsigned int));
  kb.write((const char*)&tdesc,sizeof(unsigned int));
  for(unsigned int n=0;n<nbpts;n++){//  while(true){
    KeyPoint kp(tdesc);
    kp.lireKeyLowe(kl,tdesc);
    if(kl.eof()){break;}
    kp.ecrireBin(kb);
  }
  kb.close();
  kl.close();
}

void OutilsIOKeyPoint::KeyBin2KeyLowe(const stringpp & adresse_keybin, const stringpp & adresse_keylowe){
  std::ifstream kb(adresse_keybin.c_str(),std::ios::in | std::ios::binary);
  std::ofstream kl(adresse_keylowe.c_str(),std::ios::out);
  if(!kl.good()){return;}
  unsigned int nbpts, tdesc;
  kb.read((char*)&nbpts,sizeof(unsigned int));
  kb.read((char*)&tdesc,sizeof(unsigned int));
  kl<<nbpts<<" "<<tdesc<<"\n";
  /*for(int n=0;n<nbpts;n++){//*/while(true){
    KeyPoint kp(tdesc);
    kp.lireBin(kb,tdesc);
    if(kb.eof()){break;}
    kp.ecrireKeyLowe(kl);
  }
  kb.close();
  kl.close();
}

void OutilsIOKeyPoint::KeyINRIA2KeyBin(const stringpp & adresse_keyinria, const stringpp & adresse_keybin){
  std::ifstream kl(adresse_keyinria.c_str(),std::ios::in);
  std::ofstream kb(adresse_keybin.c_str(),std::ios::out | std::ios::binary);
  if(!kl.good()){return;}
  unsigned int nbpts, tdesc;
  kl>>tdesc>>nbpts;
  kb.write((const char*)&nbpts,sizeof(unsigned int));
  kb.write((const char*)&tdesc,sizeof(unsigned int));
  for(unsigned int n=0;n<nbpts;n++){//  while(true){
    KeyPoint kp(tdesc);
    kp.lireKeyINRIA(kl,tdesc);
    if(kl.eof()){break;}
    kp.ecrireBin(kb);
  }
  kb.close();
  kl.close();
}

void OutilsIOKeyPoint::KeyBin2KeyINRIA(const stringpp & adresse_keybin, const stringpp & adresse_keyinria){
  std::ifstream kb(adresse_keybin.c_str(),std::ios::in | std::ios::binary);
  std::ofstream kl(adresse_keyinria.c_str(),std::ios::out);
  if(!kl.good()){return;}
  unsigned int nbpts, tdesc;
  kb.read((char*)&nbpts,sizeof(unsigned int));
  kb.read((char*)&tdesc,sizeof(unsigned int));
  kl<<tdesc<<" "<<nbpts<<"\n";
  /*for(int n=0;n<nbpts;n++){//*/while(true){
    KeyPoint kp(tdesc);
    kp.lireBin(kb,tdesc);
    if(kb.eof()){break;}
    kp.ecrireKeyINRIA(kl);
  }
  kb.close();
  kl.close();
}

void OutilsIOKeyPoint::KeyLowe2KeyINRIA(const stringpp & adresse_keylowe, const stringpp & adresse_keyinria){
  std::ifstream kl(adresse_keylowe.c_str(),std::ios::in);
  std::ofstream kb(adresse_keyinria.c_str(),std::ios::out);
  if(!kl.good()){return;}
  unsigned int nbpts, tdesc;
  kl>>nbpts>>tdesc;
  kb<<tdesc<<"\n"<<nbpts<<"\n";
  for(unsigned int n=0;n<nbpts;n++){//  while(true){
    KeyPoint kp(tdesc);
    kp.lireKeyLowe(kl,tdesc);
    if(kl.eof()){break;}
    kp.ecrireKeyINRIA(kb);
  }
  kb.close();
  kl.close();
}

void OutilsIOKeyPoint::KeyINRIA2KeyLowe(const stringpp & adresse_keyinria, const stringpp & adresse_keylowe){
  std::ifstream kb(adresse_keyinria.c_str(),std::ios::in);
  std::ofstream kl(adresse_keylowe.c_str(),std::ios::out);
  if(!kl.good()){return;}
  unsigned int nbpts, tdesc;
  kb>>tdesc>>nbpts;
  kl<<nbpts<<" "<<tdesc<<"\n";
  /*for(int n=0;n<nbpts;n++){//*/while(true){
    KeyPoint kp(tdesc);
    kp.lireKeyINRIA(kb,tdesc);
    if(kb.eof()){break;}
    kp.ecrireKeyLowe(kl);
  }
  kb.close();
  kl.close();
}



void OutilsIOKeyPoint::Keykoen2KeyBin(const stringpp & adresse_keykoen, const stringpp & adresse_keybin){
  std::ifstream kk(adresse_keykoen.c_str(),std::ios::in);
  std::ofstream kb(adresse_keybin.c_str(),std::ios::out | std::ios::binary);
  if(!kk.good()){return;}
  unsigned int nbpts,tdesc;
  std::string stmp,stmpp;
  double dtmp;
  kk>>stmp>>tdesc>>nbpts;
  if(stmp!="KOEN1"){return;}
  kb.write((const char*)&nbpts,sizeof(unsigned int));
  kb.write((const char*)&tdesc,sizeof(unsigned int));

  KeyPoint kp(tdesc);
  for(unsigned int n=0;n<nbpts;n++){
    getline(kk,stmp,';');
    std::istringstream iss1(stmp.c_str());
    iss1>>stmpp>>kp.x()>>kp.y()>>kp.sigma()>>kp.orientation();
    getline(kk,stmp,';');
    //    std::cout<<stmp<<std::endl;
    std::istringstream iss2(stmp.c_str());
    for(unsigned int k=0;k<tdesc;k++){
      iss2>>dtmp;
      kp[k]=(unsigned char)(dtmp/512.0*255.0+0.5);//necessaire sinon pb de lecture char par char... et en plus normalisation par 512...
    }
    //    std::cout<<std::endl;
    kp.ecrireBin(kb);
  }
  kb.close();
  kk.close();
}

void OutilsIOKeyPoint::KeyBin2Keykoen(const stringpp & adresse_keybin, const stringpp & adresse_keykoen, const bool avecdesc){
  std::ifstream kb(adresse_keybin.c_str(),std::ios::in | std::ios::binary);
  std::ofstream kk(adresse_keykoen.c_str(),std::ios::out);
  if(!kb.good()){return;}
  unsigned int nbpts,tdesc;
  kb.read((char*)&nbpts,sizeof(unsigned int));
  kb.read((char*)&tdesc,sizeof(unsigned int));
  kk<<"KOEN1"<<std::endl<<tdesc<<std::endl<<nbpts<<std::endl;
  for(unsigned int n=0;n<nbpts;n++){
    KeyPoint kp(tdesc);
    kp.lireBin(kb,tdesc);
    if(kb.eof()){break;}
    kk<<"<CIRCLE "<<kp.x()<<" "<<kp.y()<<" "<<kp.sigma()<<" "<<kp.orientation()<<" 0>;";
    if(avecdesc){
      for(unsigned int k=0;k<tdesc;k++){
	kk<<" "<<(int)((double)(kp[k])*512.0/255.0+0.5);
      }
    }
    kk<<";"<<std::endl;
  } 
  kb.close();
  kk.close();
}

//--------------------------------------------------------------
//Sélection des points d'intérêt compris dans une fenêtre
//--------------------------------------------------------------  
void OutilsKeyPoint::selection_fenetre(const std::vector<KeyPoint> & kp, const TRectangle2D<float> & fenetre, std::vector<KeyPoint> & kpsel, const bool strict){
  std::vector<KeyPoint>::const_iterator i;
  for(i=kp.begin();i!=kp.end();++i){
    if(strict){
      if((*i).x()>fenetre.p1().x() && (*i).y()>fenetre.p1().y() && (*i).x()<fenetre.p2().x() && (*i).y()<fenetre.p2().y()){kpsel.push_back((*i));}
    }else{
      if((*i).x()>=fenetre.p1().x() && (*i).y()>=fenetre.p1().y() && (*i).x()<=fenetre.p2().x() && (*i).y()<=fenetre.p2().y()){kpsel.push_back((*i));}
    }
  }
}

void OutilsKeyPoint::selection_fenetre(std::vector<KeyPoint> & kp, const TRectangle2D<float> & fenetre, bool strict){
  std::vector<KeyPoint> kpsel;
  OutilsKeyPoint::selection_fenetre(kp,fenetre,kpsel);
  kp=kpsel;
}


void OutilsKeyPoint::DescripteurMoyen(const std::vector<KeyPoint> & vkp, KeyPoint & kpm){
  KeyPoint::const_iterator ikp;
  std::vector<double>::iterator isiftkpm;
  std::vector<double> siftkpm((*(vkp.begin())).dim(),0.0);
  for(std::vector<KeyPoint>::const_iterator ivkp=vkp.begin();ivkp!=vkp.end();++ivkp){
    isiftkpm=siftkpm.begin();
    ikp=(*ivkp).begin();
    for(;isiftkpm!=siftkpm.end();++isiftkpm,++ikp){(*isiftkpm)=(*isiftkpm)+(double)(*ikp);}
  }
  for(isiftkpm=siftkpm.begin();isiftkpm!=siftkpm.end();++isiftkpm){(*isiftkpm)=(*isiftkpm)/(double)(vkp.size());}
  kpm=KeyPoint((*(vkp.begin())).dim());
  KeyPoint::iterator ikpm=kpm.begin();
  for(isiftkpm=siftkpm.begin();isiftkpm!=siftkpm.end();++isiftkpm,++ikpm){(*ikpm)=(unsigned char)(*isiftkpm);}
}

void OutilsKeyPoint::DescripteurMoyenGrille(const std::vector<KeyPoint> & vkp, int nc, int nl, double dc, double dl, std::vector<KeyPoint> & vkpm){
  std::vector<KeyPoint>::const_iterator ivkp;
  vkpm.resize(nc*nl);
  std::vector<std::vector<KeyPoint> > vkpgrille(nc*nl);
  for(ivkp=vkp.begin();ivkp!=vkp.end();++ivkp){
    vkpgrille[((int)((*ivkp).X()/dc))*nl+(int)((*ivkp).Y()/dl)].push_back(*ivkp);
  }
  for(int c=0;c<nc;c++){
    for(int l=0;l<nl;l++){
      int n=c*nl+l;
      DescripteurMoyen(vkpgrille[n],vkpm[n]);
      vkpm[n].X()=dc*((double)c+0.5);
      vkpm[n].Y()=dl*((double)l+0.5);
    }
  }
}



#ifdef _AVECMATIS_
//--------------------------------------------------------------
//Elimination des points d'intérêt hors du masque
//--------------------------------------------------------------  
void OutilsKeyPoint::FiltrageMasque(const std::vector<KeyPoint> & kp, const TTImage<bool> & masque, std::vector<KeyPoint> & kpsel){
  std::vector<KeyPoint>::const_iterator i;
  for(i=kp.begin();i!=kp.end();++i){
    if(masque((int)(*i).x(),(int)(*i).y())!=0){kpsel.push_back((*i));}
  }
}

void OutilsKeyPoint::FiltrageMasque(std::vector<KeyPoint> & kp, const TTImage<bool> & masque){
  std::vector<KeyPoint> kpsel;
  OutilsKeyPoint::FiltrageMasque(kp,masque,kpsel);
  kp=kpsel;
}
#endif







