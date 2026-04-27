#ifndef CCODEDTARGETREFINE_H
#define CCODEDTARGETREFINE_H

#endif // CCODEDTARGETREFINE_H

#include "CodedTarget.h"
#include "MMVII_PCSens.h"
#include "cMMVII_Appli.h"
#include "cCodedTargetDescribe.h"

namespace MMVII
{
    class cAppli_CodedTargetRefine : public cMMVII_Appli
    {
    public:
        //------
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
        //------ members
        std::string                         mFSpecName;
        std::unique_ptr<cFullSpecifTarget>  mFSpec;
        std::vector<cCdTDescr>              mVDescr;
        //------ methods
    };

    /*!
    * @brief The cCdTDiscr class stores all about coded target image discretisation
    */

    class cCdTDiscr
    {
    public:
        cCdTDiscr(cCdTDescr aDes, cSensorCamPC* aCam);
        //----- members
        std::string             mName;

        //----- methods
        int         visibility();

    private:
        //----- members
        cCdTDescr               mDescr; //-> CdT description from previous computation
        cSensorCamPC*           mCam;
        std::string             mNum;
        cPixBox<2>              mExtent;//-> extent of the CdT formed by predicted corners
        cIm2D<tU_INT1>          mIm;    //-> original image
        cDataIm2D<tU_INT1>*     mDIm;
        cIm2D<tU_INT1>          mCrop;  //-> croped from input image
        cDataIm2D<tU_INT1>*     mDCrop;
        cIm2D<tU_INT1>          mCano;  //-> canonical CdT
        cDataIm2D<tU_INT1>*     mDCano;
        cIm2D<tU_INT1>          mSimu;  //-> simul. CdT
        cDataIm2D<tU_INT1>*     mDSimu;
        tU_INT1                 mVisib; //-> visibility score

        //----- methods
        void        CamExtent();//-> compute CdT extent in original image and set mVisib to 1
        void        True();     //-> computes extent from description and extract target in mTrue
        void        GenSimul();    //-> simulate from target extent
        cPt3dr      CdT2Gnd(cPt2di aPix);
        cPt2dr      CdT2Cam(cPt2di aPix);
    };
}
