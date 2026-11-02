#include "MMVII_util_tpl.h"
#include "MMVII_Ptxd.h"
#include "MMVII_Geom2D.h"

namespace MMVII
{

class cDMEpipPairResult
{
    private:
        std::string mModeGeom;
        std::string mNameEpipIm1;
        std::string mNameEpipIm2;
        std::string mNameDisparity;
        std::string mNameCorrel;
        std::string mNameMasq;

    public:
        cDMEpipPairResult(std::string aModeGeom,
                        std::string aNameEpipIm1, 
                        std::string aNameEpipIm2,
                        std::string aNameDisparity,
                        std::string aNameCorrel,
                        std::string aNameMasq):
            mModeGeom(aModeGeom),
            mNameEpipIm1(aNameEpipIm1),     
            mNameEpipIm2(aNameEpipIm2),
            mNameDisparity(aNameDisparity),
            mNameCorrel(aNameCorrel),
            mNameMasq(aNameMasq)
        {}

        cDMEpipPairResult()
        {}

        void AddData(const cAuxAr2007 & anAux)
        {
            MMVII::AddData(cAuxAr2007("ModeGeom",anAux),mModeGeom);
            MMVII::AddData(cAuxAr2007("NameEpipIm1",anAux),mNameEpipIm1);
            MMVII::AddData(cAuxAr2007("NameEpipIm2",anAux),mNameEpipIm2);
            MMVII::AddData(cAuxAr2007("NameDisparity",anAux),mNameDisparity);
            MMVII::AddData(cAuxAr2007("NameCorrel",anAux),mNameCorrel);
            MMVII::AddData(cAuxAr2007("NameMasq",anAux),mNameMasq);
        }

        //ACCESSORS
        std::string ModeGeom() {return mModeGeom;}
        std::string NameEpipIm1() {return mNameEpipIm1;}    
        std::string NameEpipIm2() {return mNameEpipIm2;}
        std::string NameDisparity() {return mNameDisparity;}
        std::string NameCorrel() {return mNameCorrel;}
        std::string NameMasq() {return mNameMasq;}

};

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

inline void AddData(const  cAuxAr2007 & anAux, cDMEpipPairResult & aDMEpipPair) {
    aDMEpipPair.AddData(anAux);
}
}