#include "MMVII_ExifData.h"
#include "MMVII_Image2D.h"
#include "cMMVII_Appli.h"
#include "MMVII_Sensor.h"
#include "MMVII_Tpl_ElemStrToVal.h"
#include "MMVII_2Include_Serial_Tpl.h"

namespace MMVII
{

enum class eDispExif
{
    eSizeType=0,
    ePhgrExif,   // Xif tag having a photogrammetric effect
    eMainExif,
    eAllExif,
    eRawExif,
    eAllGDALInfo,
    eNbVals
};

template<> cE2Str<eDispExif>::tMapE2Str cE2Str<eDispExif>::mE2S
           {
                {eDispExif::eSizeType,"SizeType"},
                {eDispExif::ePhgrExif,"PhgrExif"},
                {eDispExif::eMainExif,"MainExif"},
                {eDispExif::eAllExif,"AllExif"},
                {eDispExif::eRawExif,"RawExif"},
                {eDispExif::eAllGDALInfo,"AllGDALInfo"}
           };

MACRO_INSTANTIATE_STRIO_ENUM(eDispExif,"DispExif")


class cAppli_ImageMetada : public cMMVII_Appli
{
public :
    cAppli_ImageMetada(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec,bool isBasic);
    int Exe() override;
    cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
    cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;

    std::vector<std::string>  Samples() const override;

private :
    /** Process one image , ie:
         - compute the hash-code if ! ForPrint
         - print else
         The hash-code will allow to group together images having same meta data
     */
    void MakeOneImage(const std::string & aName,bool ForPrint);


    void DispGDalInfo(const cDataFileIm2D  & aDataFileIm);
    void DispRawExif(const cDataFileIm2D  & aDataFileIm);
    void DispStdExif(const cDataFileIm2D  & aDataFileIm,const cMetaDataImage & aMDI);




    template <class Type> inline void FunctionDispExif(const std::string& aXifName,const Type & aXifVal);
    template <class TypeXif,class TypeMMVI2> inline void FunctionDispExif_MMVII
                                         (const std::string&aXifName,const TypeXif &aXifVal,
                                          const TypeMMVI2 &aValV2,const TypeMMVI2 & aDefV2);



    cPhotogrammetricProject  mPhgrProj;
    std::string    mNameIn;  ///< Input image name
    eDispExif      mDisp;
    cExifData      mCurXif;
    std::size_t    mHash;
    bool           mIterPrint;

    std::map<std::size_t,std::list<std::string>>  mMergedNames;
    std::map<std::string,std::size_t>             mName2Hash;

};


cAppli_ImageMetada:: cAppli_ImageMetada(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec,bool isBasic) :
    cMMVII_Appli (aVArgs,aSpec),
    mPhgrProj(*this),
    mDisp(eDispExif::eMainExif),
    mHash (0)
{
}


cCollecSpecArg2007 & cAppli_ImageMetada::ArgObl(cCollecSpecArg2007 & anArgObl)
{
    return
        anArgObl
        <<   Arg2007(mNameIn,"Name of input file",{{eTA2007::MPatFile,"0"},eTA2007::FileImage})
        ;
}

cCollecSpecArg2007 & cAppli_ImageMetada::ArgOpt(cCollecSpecArg2007 & anArgOpt)
{
    return
        anArgOpt
        <<   AOpt2007(mDisp,"Disp","What to display in enum values",{eTA2007::HDV})
        ;
}

std::vector<std::string>  cAppli_ImageMetada::Samples() const
{
    return {
           "MMVII ImageMetadata IMG_0335.JPG",
           "MMVII ImageMetadata .*JPG",
           "MMVII ImageMetadata .*JPG Disp=PhgrExif"
    };
}


template<typename T>
std::ostream& operator<<(std::ostream& os, std::optional<T> const& opt)
{
    return opt ? os << opt.value() : os << "<NULL>";
}


template <class Type>  void cAppli_ImageMetada::FunctionDispExif(const std::string& aNameKey,const Type & aValue)
{
    if (mIterPrint)
       StdOut() << Color::command << "  -(exif) "<< aNameKey << Color::end << ": "  << aValue << std::endl;
    hash_combine(mHash,aNameKey);
    hash_combine(mHash,aValue);
}

template <class TypeXif,class TypeMMVI2>
   void cAppli_ImageMetada::FunctionDispExif_MMVII
        (
             const std::string&aXifName,const TypeXif &aXifVal,
             const TypeMMVI2 &aValV2,const TypeMMVI2 & aDefV2
        )
{
       bool HasV2 = (aValV2 != aDefV2);
       bool HasXif = mCurXif.Valid() && aXifVal.has_value();
       bool V2EqXif = HasV2 && HasXif && (aXifVal.value()==aValV2);

       if ((!HasV2)  || V2EqXif)
       {
           if (mCurXif.Valid())
              FunctionDispExif(aXifName,aXifVal);
       }
       else
       {
          if (mIterPrint)
          {
              if (HasXif)
              {
                  StdOut() << Color::warning << "  *** MMVII, ModifiedByUserRule " << Color::end
                           << Color::command << aXifName << Color::end
                           << ": [ " <<  aXifVal   << "  ===>  " << aValV2 << " ]\n";
              }
              else
              {
                 StdOut() << Color::descr << "  *** MMVII, SetByUserRule " << Color::end
                          << Color::command << aXifName << Color::end
                          << ": " << aValV2 << "\n";
              }
          }
           hash_combine(mHash,MMVIIBin2007+":"+aXifName);
           hash_combine(mHash,aValV2);
       }
}


void cAppli_ImageMetada::DispGDalInfo(const cDataFileIm2D  & aDataFileIm)
{
    auto allMetadata = aDataFileIm.AllMetadata();
    for (const auto& aDomain : allMetadata )
    {
        if (mIterPrint)
            StdOut() << "- Domain : " << (aDomain.first.empty() ? "<NULL>" : "\"" + aDomain.first + "\"") << "\n";
        for (const auto& aMetadata : aDomain.second)
        {
            hash_combine(mHash,aMetadata);

            if (mIterPrint)
                StdOut() << "  . \"" << aMetadata << "\"\n";
        }

    }
}

void cAppli_ImageMetada::DispRawExif(const cDataFileIm2D  & aDataFileIm)
{
    auto anExifList = aDataFileIm.ExifStrings();
    if (anExifList.empty())
    {
        if (mIterPrint)
           StdOut() << "No Exif metadata" << std::endl;
    }
    else
    {
        for (const auto &s : anExifList)
        {
            if (mIterPrint)
               StdOut() << s << std::endl;
        }
    }
}

void cAppli_ImageMetada::DispStdExif(const cDataFileIm2D  & aDataFileIm,const cMetaDataImage & aMDI)
{
    mCurXif = aDataFileIm.ExifData();


   // MPD pas mal de remaniement pour :
   //  - conserver le No Exif Data ajoute par Jo quand il n'y a pas de xif
   //  - conserver quand meme dans ce cas les xif du a des modif utilisateur
   //  - je pense qu'a terme il faudrait "de-macroiser" ce code avec des template

#define DISP_EXIF(key) {FunctionDispExif(#key, mCurXif.m##key);}

// this macro disp standard xif, and user's changed value if apply
#define DISP_XIF_MMVII(key,val,def) {FunctionDispExif_MMVII(#key, mCurXif.m##key,val,def);}




   // DISP_EXIF(FocalLength_mm);
    DISP_XIF_MMVII(FocalLength_mm,aMDI.FocalMM(true),-1.0);

    //DISP_EXIF(FocalLengthIn35mmFilm_mm);
    DISP_XIF_MMVII(FocalLengthIn35mmFilm_mm,aMDI.FocalMMEqui35(true),-1.0);

    //DISP_EXIF(Model);
    DISP_XIF_MMVII(Model,aMDI.CameraName(true),std::string(""));
    /// the information on size of sensor are generally not provided by xif, so user
    /// has to indicate it in the data-base, we indicate if it was found
    if (mIterPrint)
    {

        std::string aNameMMVII = aMDI.CameraName(true);

        if (aNameMMVII!="")
        {
            const cElemCamDataBase * anElCDB = mPhgrProj.GetCamFromNameCam(aNameMMVII,true);
            if (anElCDB)
            {
                 StdOut()  << Color::descr << "  *** MMVII, Camera model is in data base :" << Color::end
                           << Color::command << " SzSensor=" << Color::end << anElCDB-> mSzSensor_Mm << " mm "
                            << Color::command << " SzPixel="   << Color::end << anElCDB->mSzPixel_Micron << " micron"
                           << "\n";
            }
            else
            {
                 StdOut()  << Color::warning << " !!!  Camera model is NOT in data base !!!" << Color::end << "\n";
            }
        }


    }
    if (aMDI.Has_InternalCalibGeomIdent())
    {
        std::string aNameCam = aMDI.InternalCalibGeomIdent();
        if (mIterPrint)
           StdOut() << Color::descr << "  *** MMVII, " << Color::command << " Cam-Ident " << Color::end << " : "<< aNameCam << "\n";
        hash_combine(mHash,aNameCam);
    }


   if (! mCurXif.Valid())
   {
       if (mIterPrint)
           StdOut() << Color::warning << "--- No Exif metadata---" << Color::end<< std::endl;
       return;
   }

   if (mDisp==eDispExif::ePhgrExif)
       return;


   DISP_EXIF(PixelXDimension);
   //aMDI.InternalCalibGeomIdent()
 //  const auto & aList = mMergedNames[aName];
//   mMergedNames
   DISP_EXIF(PixelYDimension);


   DISP_EXIF(FNumber);
   DISP_EXIF(ExposureTime_s);
   DISP_EXIF(Orientation);
   DISP_EXIF(Make);
   DISP_EXIF(LensMake);
   DISP_EXIF(LensModel);

   if (mDisp==eDispExif::eMainExif)
       return;


   DISP_EXIF(XResolution);
   DISP_EXIF(YResolution);
   DISP_EXIF(ResolutionUnit);
   DISP_EXIF(FocalPlaneXResolution);
   DISP_EXIF(FocalPlaneYResolution);
   DISP_EXIF(FocalPlaneResolutionUnit);

   DISP_EXIF(DateTime);
   DISP_EXIF(SubSecTime);
   DISP_EXIF(DateTimeOriginal);
   DISP_EXIF(SubSecTimeOriginal);
   DISP_EXIF(DateTimeDigitized);
   DISP_EXIF(SubSecTimeDigitized);


   const auto default_precision{std::cout.precision()};
   constexpr auto max_precision{std::numeric_limits<long double>::digits10};

   StdOut() << std::setprecision(max_precision);
   DISP_EXIF(DateTimeNumber_s);
   DISP_EXIF(DateTimeOriginalNumber_s);
   DISP_EXIF(DateTimeDigitizedNumber_s);

   DISP_EXIF(GPSLongitude_deg);
   DISP_EXIF(GPSLatitude_deg);
   DISP_EXIF(GPSAltitude_m);
 // ??  StdOut() << std::setprecision(default_precision);

   DISP_EXIF(GPSDateStamp);
   DISP_EXIF(GPSTimeStamp);
   DISP_EXIF(GPSTimeUTC_s);
   DISP_EXIF(GPSTimeUTC_ns);
   DISP_EXIF(ExifVersion);

#undef DISP_EXIF
#undef DISP_XIF_MMVII
   StdOut() << std::setprecision(default_precision);
}



void cAppli_ImageMetada::MakeOneImage(const std::string & aName,bool ForPrint)
{
    mHash = 0;

    //const auto default_precision{std::cout.precision()};
    //constexpr auto max_precision{std::numeric_limits<long double>::digits10};

    cMetaDataImage aMDI = mPhgrProj.GetMetaData(aName);

    auto aDataFileIm=cDataFileIm2D::Create(aName,eForceGray::No);

    if (mIterPrint)
    {

        {
            size_t aHash = mName2Hash[aName];
            auto & aList = mMergedNames[aHash];
            if (aList.empty())
                return;
            for (const auto & aNameEquiv : aList)
                StdOut() << Color::title << "------------- " << aNameEquiv <<"------------- " << Color::end <<  std::endl;
            aList.clear();
        }

        StdOut() << "  "
                 << Color::command << "Size: " << Color::end << aDataFileIm.Sz()
                 << Color::command << ", Type: " << Color::end  << ToStr(aDataFileIm.Type())
                 << Color::command << ", Channels: " << Color::end << aDataFileIm.NbChannel() << std::endl;


    }
    else
    {
        hash_combine(mHash,std::string("SizeType"));
        hash_combine(mHash,HashValue(aDataFileIm.Sz(),true));
        hash_combine(mHash,(int)aDataFileIm.Type());
        hash_combine(mHash,(int)aDataFileIm.NbChannel());
    }
    switch (mDisp)
    {
        case eDispExif::eSizeType:
        case eDispExif::eNbVals:
        {
            return;
        }

        case eDispExif::ePhgrExif:
        case eDispExif::eMainExif:
        case eDispExif::eAllExif:
        {
             mCurXif = aDataFileIm.ExifData();

             DispStdExif(aDataFileIm,aMDI);
             return;
        }

        case eDispExif::eRawExif:
        {
            DispRawExif(aDataFileIm);
            return;
        }

        case eDispExif::eAllGDALInfo:
        {
             DispGDalInfo(aDataFileIm);
             return;
        }
    }

}

int cAppli_ImageMetada::Exe()
{
    mPhgrProj.FinishInit();



  //  const auto default_precision{std::cout.precision()};
  //  constexpr auto max_precision{std::numeric_limits<long double>::digits10};

    for (const auto anIter : {false,true})
    {
        mIterPrint = anIter;
       for (const auto & aName : VectMainSet(0))
       {
           MakeOneImage(aName,true);
           if (!mIterPrint)
           {
               mMergedNames[mHash].push_back(aName);
               mName2Hash[aName]=mHash;
//               StdOut() << "NAME=" << aName << " H=" << mHash << "\n";
           }
       }
    }
    if (VectMainSet(0).size() > 1)
    {
        StdOut() << Color::title << "--Number of set grouped by meta data : " << Color::end << mMergedNames.size() << "\n";
    }
    return EXIT_SUCCESS;
}



static tMMVII_UnikPApli Alloc_ImageMetadata(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec)
{
    return tMMVII_UnikPApli(new cAppli_ImageMetada(aVArgs,aSpec,true));
}

cSpecMMVII_Appli  TheSpec_ImageMetadata
    (
        "ImageMetadata",
        Alloc_ImageMetadata,
        "Display Exif and Metadata from image file",
        {eApF::ImProc},
        {eApDT::Image},
        {eApDT::Console},
        __FILE__
        );


} // namespace MMVII
