#include "MMVII_ExifData.h"
#include "MMVII_Image2D.h"
#include "cMMVII_Appli.h"

namespace MMVII
{

class cAppli_ImageMetada : public cMMVII_Appli
{
public :
    cAppli_ImageMetada(const std::vector<std::string> & aVArgs,const cSpecMMVII_Appli & aSpec,bool isBasic);
    int Exe() override;
    cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override ;
    cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override ;

private :
    std::string mNameIn;  ///< Input image name
    int mDisp;
};


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
        <<   AOpt2007(mDisp,"Disp","0:Size & type, 1:Main Exif, 2:All Exif, 3:Raw Exif  4:all GDAL info ",{eTA2007::HDV,{eTA2007::Range,"[0,4]"}})
        ;
}



template<typename T>
std::ostream& operator<<(std::ostream& os, std::optional<T> const& opt)
{
    return opt ? os << opt.value() : os << "<NULL>";
}


int cAppli_ImageMetada::Exe()
{
    const auto default_precision{std::cout.precision()};
    constexpr auto max_precision{std::numeric_limits<long double>::digits10};

    for (const auto & aName : VectMainSet(0))
    {
        auto aDataFileIm=cDataFileIm2D::Create(aName,eForceGray::No);
        StdOut() << "####### " << aDataFileIm.Name() <<": " << std::endl;
        StdOut() << "Size: " << aDataFileIm.Sz() << ", Type: "  << ToStr(aDataFileIm.Type()) << ", Channels: " << aDataFileIm.NbChannel() << std::endl;
        switch (mDisp) {
        case 0:
            break;
        case 1:
        case 2:
        {
            cExifData anExif = aDataFileIm.ExifData();
            if (! anExif.Valid())
            {
                StdOut() << "No Exif metadata" << std::endl;
                break;
            }

#define DISP_EXIF(key) StdOut() << #key << ": " << anExif.m##key << std::endl;

            DISP_EXIF(PixelXDimension);
            DISP_EXIF(PixelYDimension);

            DISP_EXIF(FocalLength_mm);
            DISP_EXIF(FocalLengthIn35mmFilm_mm);
            DISP_EXIF(FNumber);
            DISP_EXIF(ExposureTime_s);
            DISP_EXIF(Orientation);
            DISP_EXIF(Make);
            DISP_EXIF(Model);
            DISP_EXIF(LensMake);
            DISP_EXIF(LensModel);

            if (mDisp == 2) {
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
        case 3:
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
        case 4:
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
    }
    return EXIT_SUCCESS;
}


cAppli_ImageMetada:: cAppli_ImageMetada(const std::vector<std::string> &  aVArgs,const cSpecMMVII_Appli & aSpec,bool isBasic) :
    cMMVII_Appli (aVArgs,aSpec),
    mDisp(1)
{
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
