#include "PreparationImage.hpp"


//=================================================
//classe PreparationImage 
// Quelques fonctions de traitement afin de preparer les images :
// - crop de la dalle
// - pretraitements (choix du canal, calcul du gradient, reetalement...) parametres par la classe ParametragePretraitement
// - conversion du format externe (classe ImageMemoire) vers le format interne (VLImageBuffer)
//-------------------------------------------------


//==============================================
// Conversion image format interne / format rxterne 
// i.e. vers VL::PgmBuffer, et inversement
//----------------------------------------------

void ConversionVLBufferVersImage(ImageMemoire & ima, VL::pixel_t const* im, const int width, const int height){
  ima=ImageMemoire(MMVII::cPt2di(width,height));
  MMVII::tREAL4** ptima=ima.DIm().ExtractRawData2D();
  for(int y=0;y<height;++y){
    for(int x=0;x<width;++x){
      ptima[y][x]=(float)(*im++);
    }
  }
  
}

void ConversionImageVersVLBuffer(/*const*/ ImageMemoire & ima, VL::pixel_t * im, const int width, const int height){
  MMVII::tREAL4** ptima=ima.DIm().ExtractRawData2D();
  for(int y=0;y<height;++y){
    for(int x=0;x<width;++x){
      (*im++)=ptima[y][x];//ima(x,y);
    }
  }
}


const float gradient(MMVII::tREAL4** pt, const int i, const int j, const int NbColonnes, const int NbLignes){
  if(i<=0 || j<=0 || i>=(int)NbColonnes-1 || j>=(int)NbLignes-1){return 0;}//bord a 0
  double val=pow((double)pt[j][i+1]-(double)pt[j][i-1],2.0)+pow((double)pt[j+1][i]-(double)pt[j-1][i],2.0);
  return val=sqrt(val)/2.0;
}



void PreparationImage::ChangementFormatExterneVersInterne(/*const*/ ImageMemoire & im,
							  VL::ImageBuffer & buffer, 
							  const int canal,//=0
							  const bool usegradient,//=false
							  const float maxim,//=-1
							  const float minim//=0
							  ) const{
  buffer.width=im.DIm().Sz().x();//im.GetNbColonnes();
  buffer.height=im.DIm().Sz().y();//GetNbLignes();
  buffer.data=new VL::pixel_t [ buffer.width*buffer.height ];
  VL::pixel_t * pt; 
  pt=buffer.data;//n=j*nbcolonnes+i
  float valmax=0;
  double coefreetal;
  float val;
  MMVII::tREAL4** ptima=im.DIm().ExtractRawData2D();
  for(int j=0;j<(int)buffer.height;j++){
    for(int i=0;i<(int)buffer.width;i++){
      if(!usegradient){
	val=ptima[j][i];//im.valeur(i,j,canal);//(float)ValeurPixel(im,i,j,canal);
      }else{
	//TODO
	val=gradient(ptima,i,j,buffer.width,buffer.height);//im.gradient(i,j,canal);//(float)Gradient(im,i,j,canal);
      }
      (*pt)=val;
      if(maxim<=0){
	if(val>valmax){valmax=val;}
      }
      pt++;
    }
  }
  //  std::cout<<"valeur maximale : "<<valmax<<std::endl;
  if(maxim<=0){
    coefreetal=1.0/(valmax-minim);
  }else{
    coefreetal=1.0/(maxim-minim);
  }
  pt=buffer.data;
  for(int n=0;n<buffer.width*buffer.height;n++,pt++){
    if(minim!=0){(*pt)=((*pt)-minim)*coefreetal;}else{(*pt)=(*pt)*coefreetal;}
  }
    
  //  ImageMemoire I; ConversionVLBufferVersImage(I,buffer.data,buffer.width,buffer.height); I.get().Sauver("tmp_new.tif");//test 
  //SauverImage(std::string("imtest.tif"),buffer.data,buffer.width,buffer.height);
}

void PreparationImage::ChangementFormatExterneVersInterne(/*const*/ ImageMemoire & im, VL::ImageBuffer & buffer, const ParametragePretraitement & parametrage) const{
  ChangementFormatExterneVersInterne(im,buffer,parametrage.canal(),parametrage.usegradient(),parametrage.maxim(),parametrage.minim());
}





  


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

const bool PreparationImage::GetDalle(const AccesImage & imgf, 
				      const ParametragePretraitement parametragepretraitement,
				      const int i0_grande_image, const int j0_grande_image, const int nc, const int nl, 
				      const int marge, int & xmin, int & ymin, int & xmax, int & ymax,
				      VL::ImageBuffer & buffer, 
				      const bool verbose) const{
  
  if(verbose) std::cout<<"Crop..."<<std::endl;
  //  std::cout<<imgf.GetTaille()<<" "<<imgf.GetNbCanaux()<<" "<<imgf.TypePixel()<<std::endl;
  int x1=i0_grande_image-marge; 
  int y1=j0_grande_image-marge;  
  //int x2=nc+i0_grande_image+marge;  
  //  int y2=nl+j0_grande_image+marge;
  //x1=0 x2=nbcolonnes
    
  //ImageMemoire Imagette;
  //if(verbose) std::cout<<"Coordonnees du crop dans la grande image : "<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<std::endl;
  //Imagette.Crop(imgf, x1, y1, x2, y2, 1);
  ImageMemoire Imagette(MMVII::cPt2di(nc+2*marge,nl+2*marge));
  Imagette.Read(imgf(),MMVII::cPt2di(x1,y1));
  
  if(verbose)std::cout<<"Emprise dans la grande image: ("<<i0_grande_image<<";"<<j0_grande_image<<")->("<<nc-1+i0_grande_image<<";"<<nl-1+j0_grande_image<<")"<<std::endl;
  //  Imagette.Sauver(TChaine("im-")+TChaine(i0)+TChaine("-")+TChaine(j0)+TChaine(".tif"));

  ChangementFormatExterneVersInterne(Imagette,buffer,parametragepretraitement);
  
  //Calcul des marges reelles dans lequel SIFT va s'appliquer (cf problemes de bord d'image)
  xmin=marge-(int)fabs((double)std::min(0,i0_grande_image-marge));
  ymin=marge-(int)fabs((double)std::min(0,j0_grande_image-marge));
  xmax=std::min(xmin+nc,(int)imgf.GetTaille().X())-1;
  ymax=std::min(ymin+nl,(int)imgf.GetTaille().Y())-1;
  //On ne conservera donc que les points Pt de Lowe  tels que xmin<=Pt.x<=xmax && ymin<=Pt.y<=ymax
  if(verbose) std::cout<<"Emprise \"utile\" au sein de la dalle : ("<<xmin<<" "<<ymin<<" "<<xmax<<" "<<ymax<<")"<<std::endl;
  if(verbose) std::cout<<"taille dalle :"<<Imagette.DIm().Sz().x()<<","<<Imagette.DIm().Sz().y()<<" autorisee : ("<<xmax-xmin+1<<","<<ymax-ymin+1<<")"<<std::endl;
  if(xmax<=xmin || ymax<=ymin || Imagette.DIm().Sz().x()==0 || Imagette.DIm().Sz().y()==0 ){xmin=-1; ymin=-1; xmax=-1; ymax=-1; return false;}
  
  if(verbose) std::cout<<"Fin du crop...OK"<<std::endl;
  return true;
}


const bool PreparationImage::GetDalle(const std::string & adresseimage, 
				      const ParametragePretraitement parametragepretraitement,
				      const int i0_grande_image, const int j0_grande_image, const int nc, const int nl, 
				      const int marge, int & xmin, int & ymin, int & xmax, int & ymax, 
				      VL::ImageBuffer & buffer, 
				      const bool verbose) const{
  return GetDalle(AccesImage(adresseimage),parametragepretraitement,i0_grande_image,j0_grande_image,nc,nl,marge,xmin,ymin,xmax,ymax,buffer,verbose);
}


