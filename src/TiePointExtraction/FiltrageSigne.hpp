#ifndef _FILTRAGESIGNE_HPP_
#define _FILTRAGESIGNE_HPP_

//=============================================================
//Classe FiltreSigne
// Filtrage des appariements en regardant le signe des kp
//=============================================================
class FiltreSigne{
 public:
  //Constructeur / Destructeur
  FiltreSigne(){};
  virtual ~FiltreSigne(){}

  void filtre(std::vector<KeyPoint> & pts1, std::vector<KeyPoint> & pts2) const{
    if(pts1.size()!=pts2.size()){return;}
    if(pts1.empty()){return;}
    std::vector<KeyPoint> ptstmp1,ptstmp2;
    ptstmp1.reserve(pts1.size());
    ptstmp2.reserve(pts2.size());
    std::vector<KeyPoint>::const_iterator i1,i2;
    i1=pts1.begin();
    i2=pts2.begin();
    for(;i1!=pts1.end();++i1,++i2){
      if(filtre((*i1),(*i2))){
	ptstmp1.push_back(*i1);
	ptstmp2.push_back(*i2);
      }
    }
    pts1=ptstmp1;
    pts2=ptstmp2;
  }

  const bool filtre(const KeyPoint & kp1, const KeyPoint & kp2) const{
    return (kp1.signe()==kp2.signe());
  }
  const bool operator()(const KeyPoint & kp1, const KeyPoint & kp2) const{
    return filtre(kp1,kp2);
  }
    
  void operator()(std::vector<KeyPoint> & pts1, std::vector<KeyPoint> & pts2) const{
    filtre(pts1,pts2);
  }

};


#endif
