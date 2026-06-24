#include "MMVII_ExifData.h"
#include "MMVII_Image2D.h"
#include "cMMVII_Appli.h"
#include "MMVII_Sensor.h"
#include "MMVII_Tpl_ElemStrToVal.h"

namespace MMVII
{

enum class eDispExif
{
    eSizeType=0,
    eMainExif,
    eAllExif,
    eRawExif,
    eAllGDALInfo,
    eNbVals
};

template<> cE2Str<eDispExif>::tMapE2Str cE2Str<eDispExif>::mE2S
           {
                {eDispExif::eSizeType,"SizeType"},
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

private :
    cPhotogrammetricProject  mPhgrProj;
    std::string mNameIn;  ///< Input image name
    eDispExif mDisp;
};


cAppli_ImageMetada:: cAppli_ImageMetada(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec,bool isBasic) :
    cMMVII_Appli (aVArgs,aSpec),
    mPhgrProj(*this),
    mDisp(eDispExif::eMainExif)
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
        <<   AOpt2007(mDisp,"Disp","What to display: SizeType|MainExif|AllExif|RawExif|AllGDALInfo",{eTA2007::HDV})
        ;
}



template<typename T>
std::ostream& operator<<(std::ostream& os, std::optional<T> const& opt)
{
    return opt ? os << opt.value() : os << "<NULL>";
}



int cAppli_ImageMetada::Exe()
{
    mPhgrProj.FinishInit();



    const auto default_precision{std::cout.precision()};
    constexpr auto max_precision{std::numeric_limits<long double>::digits10};

    for (const auto & aName : VectMainSet(0))
    {
        cMetaDataImage aMDI = mPhgrProj.GetMetaData(aName);

        auto aDataFileIm=cDataFileIm2D::Create(aName,eForceGray::No);
        StdOut() << "####### " << aDataFileIm.Name() <<": " << std::endl;
        StdOut() << "Size: " << aDataFileIm.Sz() << ", Type: "  << ToStr(aDataFileIm.Type()) << ", Channels: " << aDataFileIm.NbChannel() << std::endl;
        switch (mDisp) {
        case eDispExif::eSizeType:
        case eDispExif::eNbVals:
            break;
        case eDispExif::eMainExif:
        case eDispExif::eAllExif:
        {
            cExifData anExif = aDataFileIm.ExifData();
            if (! anExif.Valid())
            {
                StdOut() << "No Exif metadata" << std::endl;
                break;
            }

#define DISP_EXIF(key) StdOut() << #key << ": " << anExif.m##key << std::endl;
// this macro disp standard xif, and user's changed value if apply
#define DISP_XIF_MMVII(key,val,def) DISP_EXIF(key); \
    if ((val != def) && (anExif.m##key.value_or(def)!=val)  )  StdOut() << "  *** SetByUserRule " << val << "\n";


            DISP_EXIF(PixelXDimension);
            DISP_EXIF(PixelYDimension);

           // DISP_EXIF(FocalLength_mm);
            DISP_XIF_MMVII(FocalLength_mm,aMDI.FocalMM(true),-1);

            //DISP_EXIF(FocalLengthIn35mmFilm_mm);
            DISP_XIF_MMVII(FocalLengthIn35mmFilm_mm,aMDI.FocalMMEqui35(true),-1);

            DISP_EXIF(FNumber);
            DISP_EXIF(ExposureTime_s);
            DISP_EXIF(Orientation);
            DISP_EXIF(Make);

            //DISP_EXIF(Model);
            DISP_XIF_MMVII(Model,aMDI.CameraName(true),"");


            /// the information on size of sensor are generally not provided by xif, so user
            /// has to indicate it in the data-base, we indicate
            {
                std::optional<std::string> aModel = anExif.mModel;
                if (aModel.has_value())
                {
                    const cElemCamDataBase * anElCDB = mPhgrProj.GetCamFromNameCam(aModel.value(),true);
                    if (anElCDB)
                    {
                         StdOut()  << "  *** Camera model is in data base :"
                                   << " SzSensor=" << anElCDB-> mSzSensor_Mm << " mm "
                                   << " SzPixel="   << anElCDB->mSzPixel_Micron << " mu"
                                   << "\n";
                    }
                    else
                    {
                         StdOut()  << " !!!  Camera model is NOT in data base\n";
                    }
                }
            }
            DISP_EXIF(LensMake);
            DISP_EXIF(LensModel);

            if (mDisp == eDispExif::eAllExif) {
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

                StdOut() << std::setprecision(max_precision);
                DISP_EXIF(DateTimeNumber_s);
                DISP_EXIF(DateTimeOriginalNumber_s);
                DISP_EXIF(DateTimeDigitizedNumber_s);

                DISP_EXIF(GPSLongitude_deg);
                DISP_EXIF(GPSLatitude_deg);
                DISP_EXIF(GPSAltitude_m);
                StdOut() << std::setprecision(default_precision);

                DISP_EXIF(GPSDateStamp);
                DISP_EXIF(GPSTimeStamp);
                DISP_EXIF(GPSTimeUTC_s);
                DISP_EXIF(GPSTimeUTC_ns);

                DISP_EXIF(ExifVersion);
            }
#undef DISP_EXIF
            StdOut() << std::setprecision(default_precision);
            break;
        }
        case eDispExif::eRawExif:
        {
            auto anExifList = aDataFileIm.ExifStrings();
            if (anExifList.empty())
            {
                StdOut() << "No Exif metadata" << std::endl;
            } else {
                for (const auto &s : anExifList)
                    StdOut() << s << std::endl;
            }
            break;
        }
        case eDispExif::eAllGDALInfo:
        {
            auto allMetadata = aDataFileIm.AllMetadata();
            for (const auto& aDomain : allMetadata ) {
                StdOut() << "- Domain : " << (aDomain.first.empty() ? "<NULL>" : "\"" + aDomain.first + "\"") << "\n";
                for (const auto& aMetadata : aDomain.second) {
                    StdOut() << "  . \"" << aMetadata << "\"\n";
                }
            }
            break;
        }
        }
        StdOut() << std::endl;




//   const cElemCamDataBase *  cPhotogrammetricProject::GetCamFromNameCam(const std::string& aNameCam,bool SVP) const

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
