#include "MMVII_util_tpl.h"
#include "MMVII_Ptxd.h"

namespace MMVII
{
class cCloudRaster
{
    private:
        std::string mNameDem;
        std::string mNameCorrel;
        std::string mNameMasq;
        cPt2di mSzPix;
        cAffin2D<tREAL8> mAffineGeoRef;

    public:
        cCloudRaster(std::string aNameDem,
                    std::string aNameCorrel,
                    std::string aNameMasq,
                    cPt2di aSzPix,
                    cAffin2D<tREAL8> anAff):
            mNameDem(aNameDem),
            mNameCorrel(aNameCorrel),
            mNameMasq(aNameMasq),
            mSzPix(aSzPix),
            mAffineGeoRef(anAff)
        {
        }

        cCloudRaster()
        {
        }

        void AddData(const cAuxAr2007 & anAux)
        {
            //names
            MMVII::AddData(cAuxAr2007("NameDem",anAux) ,mNameDem);
            MMVII::AddData(cAuxAr2007("NameCorrel",anAux),mNameCorrel);
            MMVII::AddData(cAuxAr2007("NameMasq",anAux),mNameMasq);

            //Image size 
            //cPtxd<int,2> anImageSize= mBoxGlobOutPix.Sz();
            MMVII::AddData(cAuxAr2007("NbPixels",anAux),mSzPix);
            
            //Affine Transform for georeferencing
            MMVII::AddData(cAuxAr2007("AffineTransform",anAux),mAffineGeoRef);
        }


        //ACCESSORS
        std::string NameDem() {return mNameDem; }
        std::string NameCorrel() {return mNameCorrel ;}
        std::string NameMasq() { return mNameMasq ; }
        
        cPt2di SzPix() {return mSzPix ;}
        cAffin2D<tREAL8>  Transform() {return mAffineGeoRef; }

        
};

inline void AddData(const  cAuxAr2007 & anAux, cCloudRaster & aCldRaster) {
    aCldRaster.AddData(anAux);
}

};