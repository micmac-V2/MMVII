#ifndef CCODEDTARGETREFINE_H
#define CCODEDTARGETREFINE_H

#endif // CCODEDTARGETREFINE_H

#include "CodedTarget.h"
#include "MMVII_PCSens.h"
#include "cMMVII_Appli.h"
#include "cCodedTargetDescribe.h"
#include "MMVII_Interpolators.h"

namespace MMVII
{

    typedef cIm2D<tU_INT1>      tIm;
    typedef cDataIm2D<tU_INT1>  tDIm;

    class cAppli_CodedTargetRefine : public cMMVII_Appli
    {
    public:
        cAppli_CodedTargetRefine(const std::vector<std::string>& aVArgs,
                                 const cSpecMMVII_Appli& aSpec);
    private:
        //------ MMVII mandatory/usual stuff
        int Exe() override;
        cCollecSpecArg2007 & ArgObl(cCollecSpecArg2007 & anArgObl) override;
        cCollecSpecArg2007 & ArgOpt(cCollecSpecArg2007 & anArgOpt) override;
        cPhotogrammetricProject mPhProj;
        std::string mSpecImIn;
        bool mShow;
        bool mVisu;

        //------ members
        std::string                         mFSpecName; //-> full specification file name
        std::unique_ptr<cFullSpecifTarget>  mFSpec;     //-> full specification object
        std::vector<cAugCdT>                mVAugCdT;    //-> descriptions from Describe
        std::string                         mGlobImN;
        tIm                                 mIm;        //-> current image
        tDIm*                               mDIm;       //-> current image data
        cSensorCamPC*                       mCam;       //-> current camera
        cSetMesPtOf1Im                      mSetImMes;  //-> current image measurements
        tU_INT1                             mL1Lim;     //-> L1 limit to consider outliers from ransac TF computation

        //------ methods
        void        BuildDiscr(cCdTDiscr& aDis, bool& isOk);  //-> cCdTDiscr builder
        void        DiscrMapRefine(cCdTDiscr& aDis);                //-> cCdTDiscr CdT2Im mapping refiner
        void        Visu(cSetMesPtOf1Im& aSet);
        std::string NameVisu(const std::string & aIm, const std::string & aPref, const std::string aPost);
    };

    cPixBox<2>          BBox(std::vector<cPt2dr> aVPts, int aMin=0, int aMax=100000);//-> computes bounding box from a point vector
    std::vector<cPt2dr> Corners(const cPt2dr& aP0, const cPt2dr& aP1);//-> gets corners of a rectangle formed by aP0/aP1
    //cAff2D_r            Descr2Aff(const cCdTDescr& aDes, cSensorCamPC* aCam);//-> converts description to 2d affinity
}
