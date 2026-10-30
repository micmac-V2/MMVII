#ifndef _OUTILSSTRING__
#define _OUTILSSTRING__

#include <string>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <iostream>


#ifdef _SANSMATIS_
//#warning sansmatis
#include "ersatzmatis.h"
#else
//#warning avecmatis
#include "matis/tchaine.hpp"
#endif

namespace OutilsString{
  inline TChaine GetExtension(const TChaine & adresse);
  inline void SetExtension(TChaine & adresse, const TChaine & nouvelle_extension);
  inline TChaine newSetExtension(const TChaine & adresse, const TChaine & nouvelle_extension);
  inline void SetExtension(TChaine & adresse, const char *nouvelle_extension);
  inline TChaine newSetExtension(const TChaine & adresse, const char *nouvelle_extension);
  inline TChaine Extension(const TChaine & adresse);
  inline TChaine adressefromDirFile(const TChaine & dir, const TChaine & file);
  inline void DirFile(const TChaine & adresse, TChaine &dir, TChaine &file);
  inline TChaine File(const TChaine & adresse);
  inline TChaine Dir(const TChaine & adresse);
  inline TChaine CombineNom(const TChaine & adresse1, const TChaine & adresse2, const TChaine & nouvelle_extension="");
  inline void DeCombineNom(const TChaine & adresse, TChaine & adresse1, TChaine & adresse2, const TChaine & nouvelle_extension="");
  inline void SupprExtension(TChaine & adresse);
  inline TChaine newSupprExtension(const TChaine & adresse);
  inline TChaine int2string(int i);
  inline TChaine double2string(double d);
  inline TChaine double2string(double d, int precision);
  inline int string2int(const TChaine & str);
  inline double string2double(const TChaine & str);
  inline char Upcase(char a);
  inline TChaine newUpcase(const TChaine & str);
  inline void Upcase(TChaine & str);
};

inline void OutilsString::SupprExtension(TChaine & adresse){
  std::string::size_type pos=0;
  pos=adresse.find_last_of(".");
  if(pos!=adresse.npos){
    adresse.erase(pos,adresse.length());//pour ne pas garder de point
  }
}

inline TChaine OutilsString::newSupprExtension(const TChaine & adresse){
  TChaine newadresse(adresse);
  SupprExtension(newadresse);
  return newadresse;
}

inline TChaine OutilsString::GetExtension(const TChaine & adresse){
  TChaine ext(adresse);
  std::string::size_type pos=0;
  pos=ext.find_last_of(".");
  if(pos!=ext.npos){
    ext.erase(0,pos+1);//pour ne pas garder de point
  }  
  return ext;
}

inline void OutilsString::SetExtension(TChaine & adresse, const TChaine & nouvelle_extension){
  if(nouvelle_extension==TChaine("")){SupprExtension(adresse);return;}
  else if(nouvelle_extension==TChaine(".")){SupprExtension(adresse);adresse+=TChaine(".");return;}
#ifdef _SANSMATIS_
  std::string::size_type pos=0;
  pos=adresse.find_last_of(".");
  if(pos!=adresse.npos){
    adresse.erase(pos,adresse.length()-1);
  }
  //  if(nouvelle_extension.find_last_of(".")==nouvelle_extension.npos){
  if(nouvelle_extension.find(".")!=0){
    adresse=adresse+TChaine(".");
  }
  adresse=adresse+nouvelle_extension;
#else
  adresse.SetExtension(nouvelle_extension);
#endif
}

inline TChaine OutilsString::newSetExtension(const TChaine & adresse, const TChaine & nouvelle_extension){
  TChaine nouvelle_adresse(adresse);
#ifdef _SANSMATIS_
  SetExtension(nouvelle_adresse,nouvelle_extension);
#else
  nouvelle_adresse.SetExtension(nouvelle_extension);
#endif
  return nouvelle_adresse;
}

inline void OutilsString::SetExtension(TChaine & adresse, const char *nouvelle_extension){
#ifdef _SANSMATIS_
  SetExtension(adresse,TChaine(nouvelle_extension));
#else
  adresse.SetExtension(nouvelle_extension);
#endif
}

inline TChaine OutilsString::newSetExtension(const TChaine & adresse, const char *nouvelle_extension){
#ifdef _SANSMATIS_
  return newSetExtension(adresse,TChaine(nouvelle_extension));
#else
  TChaine nouvelle_adresse(adresse);
  nouvelle_adresse.SetExtension(nouvelle_extension);
  return nouvelle_adresse;
#endif
}

inline TChaine OutilsString::Extension(const TChaine & adresse){
#ifdef _SANSMATIS_
  TChaine ext(adresse);
  std::string::size_type pos=0;
  pos=ext.find_last_of(".");
  if(pos!=ext.npos){
    //    ext.erase(0,pos);//pour garder un point
    ext.erase(0,pos+1);//pour ne pas garder de point
    return ext;
  }else{
    return TChaine("");
  }
#else
  return adresse.Extension();
#endif
}

//Modifie le 16/04/08
inline TChaine OutilsString::adressefromDirFile(const TChaine & dir, const TChaine & file){
  std::string::size_type pos=0;
  pos=dir.find_last_of("/");
  if(pos!=dir.length()-1){//  if(pos==dir.npos){
    pos=dir.find_last_of("\\");
  }
  if(pos!=dir.length()-1){//  if(pos==dir.npos){
    return dir+TChaine("/")+file;
  }
  return dir+file;
}

inline void OutilsString::DirFile(const TChaine & adresse, TChaine &dir, TChaine &file){
#ifdef _SANSMATIS_
  std::string::size_type pos=0;
  pos=adresse.find_last_of("/");
  if(pos==adresse.npos){
    pos=adresse.find_last_of("\\");
  }
  if(pos==adresse.npos){
    file=adresse;
    dir=TChaine("");
  }else{
    file=adresse;
    dir=adresse;
    file.erase(0,pos+1);
    dir.erase(pos,adresse.length()-1);
  }
#else
  adresse.DirFile(dir,file);
#endif
}

inline TChaine OutilsString::File(const TChaine & adresse){
  TChaine file;
#ifdef _SANSMATIS_
  std::string::size_type pos=0;
  pos=adresse.find_last_of("/");
  if(pos==adresse.npos){
    pos=adresse.find_last_of("\\");
  }
  if(pos==adresse.npos){
    file=adresse;
  }else{
    file=adresse;
    file.erase(0,pos+1);
  }
#else
  TChaine tmp;
  adresse.DirFile(tmp,file);
#endif
  return file;
}

inline TChaine OutilsString::Dir(const TChaine & adresse){
  TChaine dir;
#ifdef _SANSMATIS_
  std::string::size_type pos=0;
  pos=adresse.find_last_of("/");
  if(pos==adresse.npos){
    pos=adresse.find_last_of("\\");
  }
  if(pos==adresse.npos){
    dir=TChaine("");
  }else{
    dir=adresse;
    dir.erase(pos,adresse.length()-1);
  }
#else
  TChaine tmp;
  adresse.DirFile(dir,tmp);
#endif
  return dir;
}


inline TChaine OutilsString::CombineNom(const TChaine & adresse1, const TChaine & adresse2, const TChaine & nouvelle_extension){
  TChaine adresse=adresse1;
  TChaine /*rep,fich1,*/fich2;
  //  DirFile(adresse1,rep,fich1);
  fich2=File(adresse2);
  SetExtension(adresse,".-.");
  if(nouvelle_extension!=""){
    SetExtension(fich2,nouvelle_extension);
  }
  return adresse+fich2;
}


inline void OutilsString::DeCombineNom(const TChaine & adresse, TChaine & adresse1, TChaine & adresse2, const TChaine & nouvelle_extension){
  TChaine adressetmp(File(adresse));
  SupprExtension(adressetmp);
  std::string::size_type pos=0;
  pos=adressetmp.find(".-.");
  if(pos==adressetmp.npos){
    adresse1="";
    adresse2="";
    return;
  }
  adresse1=adressetmp;
  adresse2=adressetmp;
  adresse1.erase(pos,adresse.length()-1);//  adresse1.erase(pos+1,adresse.length()-1);//ancienne version : on ajoutait un . a la fin du nom
  adresse2.erase(0,pos+3);
  //  adresse2=adresse2+".";//ancienne version : on ajoutait un . a la fin du nom
  if(nouvelle_extension!=TChaine("")){
    SetExtension(adresse1,nouvelle_extension);
    SetExtension(adresse2,nouvelle_extension);
  }
}



inline TChaine OutilsString::int2string(int i){
  // Conversion de int -> string
  std::string str;
  std::ostringstream ostr;
  ostr << i;
  str = ostr.str();
  return str;
}

inline int OutilsString::string2int(const TChaine & str){
  return atoi(str.c_str());
}

inline TChaine OutilsString::double2string(double d){
  // Conversion de int -> string
  std::string str;
  std::stringstream sstr;
  sstr << d;
  str = sstr.str();
  return str;
}

inline TChaine OutilsString::double2string(double d, int precision){
  // Conversion de int -> string
  std::string str;
  std::stringstream sstr;
  sstr <<std::fixed<<std::setprecision(precision)<< d;
  str = sstr.str();
  return str;
}


inline double OutilsString::string2double(const TChaine & str){
  return atof(str.c_str());
}

inline char OutilsString::Upcase(char a){
  if (a>='a' && a<='z') return ('A'+(a-'a')) ;
  return a ;
}

inline TChaine OutilsString::newUpcase(const TChaine & str){
  std::string strtmp(str);
  for(unsigned int l=0; l<str.length(); l++){
		strtmp[l]=Upcase (strtmp[l]);
  }
  return TChaine(strtmp);
}

inline void OutilsString::Upcase(TChaine & str){
	str=OutilsString::newUpcase(str);
}

//Suppression des char c avant et après...
inline const std::string supprimeavantapres(const std::string & s, const char asupprimer){
  if(s.empty()){return s;}
  std::string::size_type posd,posf;
  posd=s.find_first_not_of(asupprimer);if(posd==std::string::npos){posd=0;}
  posf=s.find_last_not_of(asupprimer);if(posf==std::string::npos){posf=s.length()-1;}
  return s.substr(posd,(int)posf-(int)posd+1);
}
//Suppression des elements contenus dans c avant et après...
inline const std::string supprimeavantapres(const std::string & s, const std::string & asupprimer){
  if(s.empty()){return s;}
  std::string::size_type posd,posf;
  posd=s.find_first_not_of(asupprimer);if(posd==std::string::npos){posd=0;}
  posf=s.find_last_not_of(asupprimer);if(posf==std::string::npos){posf=s.length()-1;}
  return s.substr(posd,(int)posf-(int)posd+1);
}
//Suppression des espaces ' ','\n','\t' avant et après...
inline const std::string supprimeavantapres(const std::string & s){
  if(s.empty()){return s;}
  std::string::size_type posd,posf;
  std::string asupprimer(" \n\t");
  posd=s.find_first_not_of(asupprimer);if(posd==std::string::npos){posd=0;}
  posf=s.find_last_not_of(asupprimer);if(posf==std::string::npos){posf=s.length()-1;}
  return s.substr(posd,(int)posf-(int)posd+1);
}

#ifndef __WINDOWS__
inline const std::string getnextbalisexml(std::ifstream & f){
  int lmax=200;
  char tmp[lmax];
  f.getline(tmp,lmax,'>');
  TChaine s(tmp);//std::string s(tmp);
  s=supprimeavantapres(s);
  OutilsString::Upcase(s);
  return s;
}
inline const std::string getnextentrebalisexml(std::ifstream & f){
  int lmax=200;
  char tmp[lmax];
  f.getline(tmp,lmax,'<');
  std::string s(tmp);
  s=supprimeavantapres(s);
  return s;
}
#endif

#endif
