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
    tU_INT1 MaskOutV = 255, MaskInV = 0;

    typedef cIm2D<tU_INT1>      tIm;
    typedef cDataIm2D<tU_INT1>  tDIm;

    class cCdTDiscr;
    struct cRansacSol;

    cPixBox<2> BBox(std::vector<cPt2dr> aVPts, int aMin=0, int aMax=100000);
    cRansacSol RansacLTF(std::vector<cPt2dr> aVBPts, std::vector<cPt2dr> aVWPts, tIm& aIm1, tIm& aIm2,
                        tDIm* aDMask=nullptr, tU_INT1 aMaskOutV=255, int aIt=200, int aRDist=50);
    std::vector<cPt2dr> Corners(const cPt2dr& aP0, const cPt2dr& aP1);

    cAff2D_r Descr2Aff(const cCdTDescr& aDes, cSensorCamPC* aCam);

    struct cRansacSol
    {
        cRansacSol(cPt2dr aSol, tREAL8 aL1Score);
        cPt2dr mSol;
        tREAL8 mL1Score;
    };

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
        bool mVisu;
        //------ members
        std::string                         mFSpecName;
        std::unique_ptr<cFullSpecifTarget>  mFSpec;
        std::vector<cCdTDescr>              mVDescr;
        tIm                                 mIm;        //-> current image
        tDIm*                               mDIm;       //-> current image data
        cSensorCamPC*                       mCam;       //-> current camera
        cSetMesPtOf1Im                      mSetImMes;  //-> current image measurements
        //------ methods
        std::string nameVisu();
        void        BuildDiscr(cCdTDiscr& aDis, cAff2D_r aCdT2Im);
        void        DiscrMapRefine(cCdTDiscr& aDis);
    };

    /*!
    * @brief The cCdTDiscr class stores all about coded target image discretisation
    */

    class cCdTDiscr
    {
    public:

        cCdTDiscr(const std::string& aName, const std::string& aImName);
        //----- members
        const std::string             mName;
        const std::string             mImName;

        //----- methods
        tREAL8                      RansacTFOnBits();
        tU_INT1                     CornersOnIm();
        void                        SaveCrop(const std::string& aDir);
        void                        SaveMask(const std::string& aDir);
        void                        SaveSample(const std::string& aDir);
        void                        SetCdT2Im(cAff2D_r aCdT2Im);
        cPt2dr                      CdT2Im(cPt2dr aPt, bool inverse=false);
        std::vector<cPt2dr>         VCdT2Im(std::vector<cPt2dr> aVPts, bool inverse=false);
        void                        SetExtent(cRect2 aExt);
        cRect2                      Extent();
        void                        SetMask(tIm& aMask);
        void                        SetCdT(tIm& aCdT);
        void                        SetCrop(tIm& aCrop);
        void                        Sample();
        tIm&                        Mask();
        tIm&                        Crop();
        tIm&                        CdT();
        tIm&                        Samp();
        void                        MapRefine();

    private:
        //----- members
        std::string             mNum;
        cPixBox<2>              mExtent;//-> extent of the CdT formed by predicted corners
        std::vector<cPt2dr>     mVBitCenters;
        tIm                     mIm;    //-> original image
        tDIm*                   mDIm;
        tIm                     mCrop;  //-> croped from input image
        tDIm*                   mDCrop;
        tIm                     mCdT;
        tIm                     mSamp;
        tIm                     mMask;  //-> simul. CdT
        tDIm*                   mDMask;
        tU_INT1                 mVisib; //-> visibility score
        std::vector<cPt2dr>     mVImCorners;
        cAff2D_r                mCdT2Im;
        cAff2D_r                mIm2CdT;
        std::vector<cPt2di>     mVCorners;

        //----- methods
        void            SaveIm(tDIm* aDIm, std::string aPath);
    };

}
