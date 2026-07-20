#ifndef CCODEDTARGETREFINE_H
#define CCODEDTARGETREFINE_H

#endif // CCODEDTARGETREFINE_H

#include "CodedTarget.h"
#include "MMVII_PCSens.h"
#include "cMMVII_Appli.h"
#include "cCodedTargetDescribe.h"

namespace MMVII
{

    typedef cIm2D<tU_INT1>      tIm;
    typedef cDataIm2D<tU_INT1>  tDIm;
    typedef cPixBox<2> tRect2;

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
        cSetOfAugCdt                        mAugSet;    //-> augmented targets set
        std::string                         mGlobImN;   //-> global image name
        tIm                                 mIm;        //-> current image
        tDIm*                               mDIm;       //-> current image data
        cSensorCamPC*                       mCam;       //-> current camera
        cSetMesPtOf1Im                      mSetImMes;  //-> current image measurements
        tU_INT1                             mL1Lim;     //-> L1 limit to consider outliers from ransac TF computation
        int                                 mMaskDil;   //-> inlier mask dilatation (wrt Ref image)
        std::string mRefine;//-> refine method to choose (corr, lsm)
        bool mMissedOnly;//-> only refine missed ctd

        //------ methods
        void        BuildDiscr(cCdTDiscr& aDis, bool& isOk);  //-> cCdTDiscr builder
        void        DiscrMapRefine(cCdTDiscr& aDis);                //-> cCdTDiscr CdT2Im mapping refiner
        cPt2dr      CorrelCropSamp(cCdTDiscr& aDis, bool &isOk);
        void        Visu(cSetMesPtOf1Im& aSet);
        std::string NameVisu(const std::string & aIm, const std::string & aPref, const std::string aPost);
        };

    tRect2          BBox(std::vector<cPt2dr> aVPts, int aMin=0, int aMax=100000);//-> computes bounding box from a point vector
    std::vector<cPt2dr> Corners(const cPt2dr& aP0, const cPt2dr& aP1);//-> gets corners of a rectangle formed by aP0/aP1


    class cPixBBox : public cPixBox<2>
    {
    public:
        typedef cPtxd<tINT4,2> tPt;
        cPixBBox(std::vector<tPt> aVPts);
        cPixBBox(cPixBox<2>);
    };

    template <class Type> class cOptCorrelThIm : public cDataMapping<tREAL8, 2, 1>
    {
    public:
        typedef  cDataIm2D<Type> tDIm;
        cOptCorrelThIm(tDIm& aTheorDIm, tDIm& aGlobDIm, cDataIm2D<tU_INT1>& aMaskDIm, cPixBBox aBBox);
        cPt1dr Value(const cPt2dr& aPt) const override;//-> correl score from aPt position
    private:
        tDIm& mThDIm;//-> theoretical image of the CdT generated from detected deformation
        tDIm& mGDIm;//-> global image
        tDIm& mDMask;//-> mask for correlation computation
        cPixBBox mBBox;//-> bbox of predicted cdt wrt global image
        cPt2dr mP0;//-> initial bbox up corner to compute translation
    };


    template <class Type> class cCdtIm
    {
    public:
        typedef cIm2D<Type> tIm;
        typedef cDataIm2D<Type> tDIm;
        cCdtIm();
    private:
        std::string mName;
        cPt2dr mC;//-> center of the cdt wrt global image
        tIm mIm;//-> cdt image
        tDIm mDIm;
    };

    class cCdtSampler
    {
    public:
        typedef cAffin2D<tREAL8> tAff2d;
        cCdtSampler(std::shared_ptr<cFullSpecifTarget> aFSpec);
        tIm Sample(std::string aName, tAff2d aRef2GIm, cPixBBox aBox, tIm& aRes, tIm& aMask);
    private:
        std::shared_ptr<cFullSpecifTarget> mFSpec;
    };

    //cAff2D_r            Descr2Aff(const cCdTDescr& aDes, cSensorCamPC* aCam);//-> converts description to 2d affinity
}
