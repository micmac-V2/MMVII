#include "MMVII_Sensor.h"
#include "MMVII_PCSens.h"
#include "MMVII_2Include_Serial_Tpl.h"

namespace MMVII
{

/*
  Vexcel example 4 first test
      30975 pixel ~ 124 mm
      Sz =  124/30975  =  4.003228 ~ 4  micron
      NbPix = 26460 17004
      SzMm  = 105.840 x 68.016
      Focal = 123.9
*/

/*  ******************************************************* */
/*                                                          */
/*                     cElemCamDataBase                     */
/*                                                          */
/*  ******************************************************* */

void cElemCamDataBase::AddData(const cAuxAr2007 & anAux)
{
    MMVII::AddData(cAuxAr2007("Name",anAux),mName);
    MMVII::AddData(cAuxAr2007("SzPix_micron",anAux),mSzPixel_Micron);
    MMVII::AddData(cAuxAr2007("SzSensor_mm",anAux),mSzSensor_Mm);
    MMVII::AddData(cAuxAr2007("NbPixels",anAux),mNbPixels);
}
void AddData(const cAuxAr2007 & anAux,cElemCamDataBase & anElem)
{
    anElem.AddData(anAux);
}

/*  ******************************************************* */
/*                                                          */
/*                     cCamDataBase                         */
/*                                                          */
/*  ******************************************************* */

void cCamDataBase::AddData(const cAuxAr2007 & anAux)
{
    MMVII::AddData(cAuxAr2007("CamDataBase",anAux),mMap);
}

void AddData(const cAuxAr2007 & anAux,cCamDataBase & aBase)
{
    aBase.AddData(anAux);
}

const std::map<std::string,cElemCamDataBase>  &  cCamDataBase::Map() const {return mMap;}
std::map<std::string,cElemCamDataBase>  &        cCamDataBase::Map()       {return mMap;}

/*  ******************************************************* */
/*                                                          */
/*                     cPhotogrammetricProject              */
/*                                                          */
/*  ******************************************************* */
static const std::string TheNameDBCam = "CameraDataBase.xml";

bool  cPhotogrammetricProject::OneTestMakeCamDataBase(const std::string & aDir,cCamDataBase & aDB)
{
  std::string aName = aDir + TheNameDBCam;
  if (! ExistFile(aName))
  {
      return false;
  }

  cCamDataBase aDBNew;
  ReadFromFile(aDBNew,aName);

  for (const auto & [aName,aCam] : aDBNew.Map())
  {
        aDB.Map()[aName] = aCam;
  }

  return true;
}

void cPhotogrammetricProject::MakeCamDataBase()
{
    for (int aNumM=(int)eTypeDBCam::eNbVals-1 ; aNumM>=0 ; aNumM--)
    {
        eTypeDBCam  aMode = (eTypeDBCam)  aNumM;
        std::string aDir = DirCamDataBase(aMode);

        OneTestMakeCamDataBase(aDir,mCamDataBase);
    }
}

const cElemCamDataBase *  cPhotogrammetricProject::GetCamFromNameCam(const std::string& aNameCam,bool SVP) const
{
    auto anIter = mCamDataBase.Map().find(aNameCam);
    if (anIter == mCamDataBase.Map().end())
    {
        if (! SVP)
        {
            MMVII_UnclasseUsEr("Cannot get camera in data base for " + aNameCam);
        }
        return nullptr;
    }

    return &(anIter->second);
}

std::string cPhotogrammetricProject::FileCamDataBase(eTypeDBCam aType)
{
    return DirCamDataBase(aType) + TheNameDBCam ;
}


std::string cPhotogrammetricProject::DirCamDataBase(eTypeDBCam aType)
{
    //std::string aNameLoc = "CameraDataBase.xml";

     mDPMetaData.SetDirInIfNoInit("Std");
    //  cMMVII_Appli::DirRessourcesMMVII();
     // std::string aDirLoc = mDPMetaData.DirIn();

     switch (aType)
     {
        case  eTypeDBCam::eLocalFolder :
            return mDPMetaData.FullDirIn()  ;

       case  eTypeDBCam::eLocalUser :
         return mAppli.DirUserProfile();

       case  eTypeDBCam::eLocalMVVI :
           return cMMVII_Appli::DirLocalParameters() ;

       case eTypeDBCam::eGlobalMMVII :
            return   cMMVII_Appli::DirRessourcesMMVII() ;

        default :
             MMVII_INTERNAL_ERROR("Value not handled in DirCamDataBase");
     }
     return "";
}


cPerspCamIntrCalib * cPhotogrammetricProject::GetCalibInit
                     (
                          const std::string& aNameIm,
                          eProjPC aProj,
                          const cPt3di & aDeg,
                          cPt2dr  aPP,
                          bool SVP,
                          bool isFraserModel
                     )
{
    thread_local static std::map<std::string,cPerspCamIntrCalib *> TheMapRes;

    // extract metadata : focal+name of camera
    cMetaDataImage  aMTD = GetMetaData(aNameIm);
    // if already exist return it
    std::string  anIdent = aMTD.InternalCalibGeomIdent();

    // StdOut() << "anIdentanIdent " << anIdent << " for " << aNameIm << "\n";


    cPerspCamIntrCalib  * & aRes = TheMapRes[anIdent];
    if (aRes!=nullptr)
       return aRes;

    // extract Camera from Data Base
    const cElemCamDataBase * aCam = GetCamFromNameCam(aMTD.CameraName(),SVP);
    if (aCam==nullptr)
    {
       return nullptr;
    }

    // [1] extract number of pixel
    cPt2di aNbPix  =   aMTD.NbPixels(true);  // if number of pixel was set by user or is in meta-data it has the prioriy
    cPt2di aNbPixCam(-1,-1);  // if number of pixel was set by user or is in meta-data it has the prioriy
    if ((aCam!=nullptr) && (aCam->mNbPixels.x()>0))
        aNbPixCam = aCam->mNbPixels;

    if (aNbPix.x() <=0)
    {
        // if the file exist  read the pixel from image
        if (ExistFile(aNameIm))
        {
            cDataFileIm2D aDF2 = cDataFileIm2D::Create(aNameIm,eForceGray::No);
            aNbPix =  aDF2.Sz();
        }
        else
        {
            // else try to get the number of pixel from the camera
            if (aNbPixCam.x() <=0)
                MMVII_UnclasseUsEr("Cannot compute number of pixel for " + aNameIm);
             aNbPix = aNbPixCam;
        }
    }

    tREAL8 aDownScaleIm = Norm2(aNbPixCam) / Norm2(aNbPix);

    // [2] extract the focal length in pixel
    tREAL8 aFocPix = aMTD.FocalPixel(true); // if it was explicitely set, it's the rule that applies
    if (aFocPix<0)
    {
        tREAL8 aFoc35 = aMTD.FocalMMEqui35(true);

      //  StdOut()  << "FOCEQUI3555555 " << aFoc35 << "\n";
        if (aFoc35>0)
        {
           // in fact I am not sure of foc equi 35, but I think it is suited for a 24x36 ..
           aFocPix =  NormInf(aNbPix) * (aFoc35 / 35.0);
        }
        else
        {
           if ((aCam==nullptr) || (aCam->mSzPixel_Micron.x()<=0))
              MMVII_UnclasseUsEr("Cannot compute focal in pixel for : " + aNameIm);
           cPt2dr aSzPMicron = aCam->mSzPixel_Micron;
           if (aSzPMicron.x() != aSzPMicron.y())
              MMVII_UnclasseUsEr("Non squared pixel non handled for now");
    
           aFocPix = (aMTD.FocalMM() / (aDownScaleIm*aSzPMicron.x())) * 1000.0;
        }
    }

    // If a PP in pixel was set, it must be used
    bool PPIsRel = true;
    cPt2dr aPPPix =  aMTD.PPPixel(true); // if it was explicitely set, it's the rule that applies
    if (aPPPix.x() > 0)
    {
       PPIsRel = false;
       aPP = aPPPix;
    }

    cDataPerspCamIntrCalib  aDataPCIC(anIdent, aProj, aDeg,aFocPix,aNbPix,PPIsRel,aPP,-1,isFraserModel);
    aRes = new cPerspCamIntrCalib(aDataPCIC);

    cMMVII_Appli::AddObj2DelAtEnd(aRes);

    return aRes;
}

/* ************************************************************* */
/*                                                               */
/*                   cAppli_CreateCalib                          */
/*                                                               */
/* ************************************************************* */

class cAppli_CreateCalib : public cMMVII_Appli
{
     public :
        cAppli_CreateCalib(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec);
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;
        // std::vector<std::string>  Samples() const override;
     private :

        cPhotogrammetricProject   mPhProj;
        std::string               mSpecIm;
        eProjPC                   mProj;
        cPt3di                    mDegree;
        cPt2dr                    mPPRel;
        bool                      mSystCyl;
};


cAppli_CreateCalib::cAppli_CreateCalib(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec) :
   cMMVII_Appli  (aVArgs,aSpec),
   mPhProj       (*this),
   mProj         (eProjPC::eStenope),
   mDegree       (3,1,1),
   mPPRel        (0.5,0.5),
   mSystCyl      (false)
{
}

cCollecSpecArg2007 & cAppli_CreateCalib::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return anArgObl
              <<  Arg2007(mSpecIm ,"Image filename",{{eTA2007::MPatFile,"0"},{eTA2007::FileDirProj}})
              <<  mPhProj.DPOrient().ArgDirOutMand()
           ;
}

cCollecSpecArg2007 & cAppli_CreateCalib::ArgOpt(cCollecSpecArg2007 & anArgObl)
{
  return      anArgObl
            << AOpt2007(mProj,"Proj","Projection mode ",{eTA2007::HDV})
            << AOpt2007(mDegree,"Degree","Degree for distorsion param",{{eTA2007::HDV}})
            << AOpt2007(mSystCyl,"SystCyl","Use SIA/SytCyl instead of Fraser Model",{{eTA2007::HDV},{eTA2007::Tuning}})
            // << AOpt2007(mNameBloc,"NameBloc","Set the name of the bloc ",{{eTA2007::HDV}})
    ;

}

int cAppli_CreateCalib::Exe()
{
    mPhProj.FinishInit();

    for (const auto & aNameIm : VectMainSet(0))
    {
        cPerspCamIntrCalib * aCalib = mPhProj.GetCalibInit(aNameIm,mProj,mDegree,mPPRel,false,!mSystCyl);
        mPhProj.SaveCalibPC(*aCalib);
    }

    return EXIT_SUCCESS;
}

tMMVII_UnikPApli Alloc_CreateCalib(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppli_CreateCalib(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_CreateCalib
(
     "OriCreateCalib",
      Alloc_CreateCalib,
      "Create initial internal calibration",
      {eApF::SysCo,eApF::Ori},
      {eApDT::Ori},
      {eApDT::Ori},
      __FILE__
);


/* ************************************************************* */
/*                                                               */
/*                   cAppli_AddCamInDataBas                      */
/*                                                               */
/* ************************************************************* */

class cAppli_AddCamInDataBase : public cMMVII_Appli
{
     public :
        cAppli_AddCamInDataBase(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec);

        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;
        std::vector<std::string>  Samples() const override;
     private :

        cPhotogrammetricProject   mPhgrProj;

        std::string               mNameCamera;
        eTypeDBCam                mTypeDB;
        int                       mMode;
        std::vector<tREAL8>       mParamCam;


};

cAppli_AddCamInDataBase::cAppli_AddCamInDataBase(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec):
    cMMVII_Appli(aVArgs,aSpec),
    mPhgrProj(*this)
{
}

cCollecSpecArg2007 & cAppli_AddCamInDataBase::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return anArgObl
              <<  Arg2007(mNameCamera ,"Name of Camera Model, or name image containg metadata")
              <<  Arg2007(mTypeDB ,"Which  database (in enum value)")
              <<  Arg2007(mMode ,"0:Info , 1:AddIfNew ; 2:Overwrite ; -1:supress")

           ;
}

cCollecSpecArg2007 & cAppli_AddCamInDataBase::ArgOpt(cCollecSpecArg2007 & anArgObl)
{
  return      anArgObl
            << AOpt2007(mParamCam,"Param","SzPix, SzSens, SzIm [SzPix.x,SzPix.y,SzSens.x,SzSens.y,SzIm.x,SzIm.y]",{{eTA2007::ISizeV,"[6,6]"}})
    ;

}

int cAppli_AddCamInDataBase::Exe()
{

    mPhgrProj.FinishInit();

    if (IsNameFileImage(mNameCamera))
    {
       cExifData anExif =   cExifData::CreateFromFile(mNameCamera);
       auto aModCam = anExif.mModel;
       if (aModCam.has_value())
       {
           mNameCamera = aModCam.value();
           StdOut() << "\n";
           StdOut() << "  *  Found  Exif Name Camera, use it " << mNameCamera << "\n";

           auto aF35 = anExif.mFocalLengthIn35mmFilm_mm;
           auto aFoc = anExif.mFocalLength_mm;
           auto aNbPix_X = anExif.mPixelXDimension;
           auto aNbPix_Y = anExif.mPixelYDimension;

           if (aF35.has_value() && aFoc.has_value() && aNbPix_X.has_value() && aNbPix_Y.has_value())
           {
               cPt2di aNbPix(aNbPix_X.value(),aNbPix_Y.value());
               tREAL8 aRatio = aFoc.value() / aF35.value();
               cPt2dr aSzSensor = cPt2dr(36,24) * aRatio;

               tREAL8 aSzPix = std::sqrt(MulCoord(aSzSensor) / MulCoord(aNbPix)  );

               aSzSensor  = ToR(aNbPix) * aSzPix;
               if (!IsInit(&mParamCam))
               {
                  std::vector<tREAL8> aVCalc
                                   {
                                       aSzPix*1e3,aSzPix*1e3,
                                       aSzSensor.x(),aSzSensor.y(),
                                       tREAL8(aNbPix.x()),tREAL8(aNbPix.y())
                                   };
                  SetIfNotInit(mParamCam,aVCalc);

                   StdOut() << "  *  Found exif senssor size "<<  aVCalc << "\n\n";
               }
           }
       }
      //  cMetaDataImage aMDI =  mPhgrProj.GetMetaData(mNameCamera);
    }

    std::string  aNameFile = mPhgrProj.FileCamDataBase(mTypeDB) ;
    cCamDataBase aCamDB;
    cElemCamDataBase aNewElem;

    if (mMode==0)
    {
        const cElemCamDataBase *  anElem = mPhgrProj.GetCamFromNameCam(mNameCamera,true);
        if (anElem==nullptr)
        {
            StdOut() << mNameCamera << " is NOT in any data base\n";
        }
        else
        {
            StdOut() << mNameCamera << " is at least  in one of the data base, with following information \n";
            StdOut()  << "  * PixSz="  << anElem->mSzPixel_Micron << " microns\n";
            StdOut()  << "  * SensoSz="  << anElem->mSzSensor_Mm << " mm\n";
            StdOut()  << "  * NbPixels="  << anElem->mNbPixels << " pixels\n";

        }
    }

    if (mMode>0)
    {
        if (!IsInit(&mParamCam))
        {
             MMVII_UnclasseUsEr("ParamCam must be initialized in write modes");
        }
        aNewElem.mName = mNameCamera;
        aNewElem.mSzPixel_Micron = cPt2dr(mParamCam.at(0),mParamCam.at(1));
        aNewElem.mSzSensor_Mm    = cPt2dr(mParamCam.at(2),mParamCam.at(3));
        aNewElem.mNbPixels       = cPt2di(mParamCam.at(4),mParamCam.at(5));
    }


    if (ExistFile(aNameFile))
    {
        ReadFromFile(aCamDB,aNameFile);
    }


    StdOut() <<  " FILE FOR DATA BASE =  " << mPhgrProj.FileCamDataBase(mTypeDB) << "\n";

    auto  anIter = aCamDB.Map().find(mNameCamera);

    if (anIter==aCamDB.Map().end())
    {
        StdOut() << mNameCamera << " is NOT in the data base " << E2Str(mTypeDB) << "\n";
        if (mMode>0)
        {
            aCamDB.Map()[mNameCamera] = aNewElem;
        }

    }
    else
    {
         StdOut() << mNameCamera << " is in the data base  " << E2Str(mTypeDB) << "\n";
         if (mMode==2)
         {
             aCamDB.Map()[mNameCamera] = aNewElem;
         }
         else if (mMode==-1)
         {
             aCamDB.Map().erase(anIter);
         }
    }

    if (mMode==0)
        return EXIT_SUCCESS;

    SaveInFile(aCamDB,aNameFile);


    return EXIT_SUCCESS;
}

std::vector<std::string>  cAppli_AddCamInDataBase::Samples() const
{
    return
    {
        "MMVII EditCamDataBase \"Canon EOS 5D Mark II\" LocalUser 1 Param=[6,6,36,24,5616,3744]"
    };
}

tMMVII_UnikPApli Alloc_AddCamInDataBase(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec)
{
   return tMMVII_UnikPApli(new cAppli_AddCamInDataBase(aVArgs,aSpec));
}

cSpecMMVII_Appli  TheSpec_AddCamInDataBase
(
     "EditCamDataBase",
      Alloc_AddCamInDataBase,
      "Create initial internal calibration",
      {eApF::Ori},
      {eApDT::ToDef},
      {eApDT::Xml},
      __FILE__
);

};


