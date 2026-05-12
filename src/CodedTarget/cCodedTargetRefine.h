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
    tU_INT1 MaskOutV = 255, MaskInV = 0;//-> Val(aPix) = MaskOutV i.e aPix is out of the mask area

    typedef cIm2D<tU_INT1>      tIm;
    typedef cDataIm2D<tU_INT1>  tDIm;

    class cCdTDiscr;
    struct cRansacSol;

    /*!
     * @brief RansacATF : computes Affine Transfert Function between 2 images
     * @param aVBPts    : subset of classified points as black values
     * @param aVWPts    : subset of classified points as white values
     * @param aDIm1     : input image
     * @param aDIm2     : output image
     * @param aDMask    : mask image (no residual computation on points in the mask)
     * @param aMaskInV  : value of pixels in the mask
     * @param aIt       : number of iterations
     * @param aRDist    : minimum grey level distance between two w/b subset points
     * @return
     */

    cRansacSol RansacATF(std::vector<cPt2dr> aVBPts, std::vector<cPt2dr> aVWPts, tDIm* aDIm1, tDIm* aDIm2,
                        tDIm* aDMask=nullptr, tU_INT1 aGVT=20, tU_INT1 aMaskOutV=255, int aIt=200, int aRDist=50);


    /*!
     * @brief The cRansacSol class : storage of RansacATF solution
     */

    struct cRansacSol
    {
        cRansacSol();
        cRansacSol(cPt2dr aSol);
        cRansacSol(cPt2dr aSol, cStdStatRes aL1Score, tIm aInlMask);
        cPt2dr      mSol;//-> best fit solution
        cStdStatRes mL1Score;//-> vector of inliers residuals
        tIm         mInlMask;//-> inliers mask
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
        std::string                         mFSpecName; //-> full specification file name
        std::unique_ptr<cFullSpecifTarget>  mFSpec;     //-> full specification object
        std::vector<cCdTDescr>              mVDescr;    //-> descriptions from Describe
        tIm                                 mIm;        //-> current image
        tDIm*                               mDIm;       //-> current image data
        cSensorCamPC*                       mCam;       //-> current camera
        cSetMesPtOf1Im                      mSetImMes;  //-> current image measurements
        tU_INT1                             mL1Lim;     //-> L1 limit to consider outliers from ransac TF computation

        //------ methods
        void        BuildDiscr(cCdTDiscr& aDis, cAff2D_r aCdT2Im);  //-> cCdTDiscr builder
        void        DiscrMapRefine(cCdTDiscr& aDis);                //-> cCdTDiscr CdT2Im mapping refiner
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
        bool                        ReqCt(std::string aRCt);//-> request content as a prerequisites check

        void                        Sample();//-> creates mSamp as a sampling of mCdT within mExtent
        void                        RansacTFSm2Cr();//-> computes mTFSm2Cr as a transfert function from mSamp to mCrop
        void                        VisuRansacTFSm2Cr(const std::string& aDir);//-> visu version

        cPt2dr                      Ref2Im(cPt2dr aPt, bool inverse=false);
        std::vector<cPt2dr>         VRef2Im(std::vector<cPt2dr> aVPts, bool inverse=false);

        void                        SetRef2Im(cAff2D_r aRef2Im);
        void                        SetExtent(cRect2 aExt);
        void                        SetMask(tIm& aMask);
        void                        SetRef(tIm& aRef);
        void                        SetCrop(tIm& aCrop);
        void                        SetWBCenters(const std::vector<cPt2dr>& aV);

        cRect2                      Extent();
        tIm&                        Mask();
        tIm&                        Crop();
        tIm&                        Ref();
        tIm&                        Samp();
        std::vector<cPt2dr>         SampC(char aC);//-> returns w/b bit centers w.r.t samp image coordinates
        cRansacSol                  TFSm2Cr();
        tIm                         ResTFSm2Cr();

        void                        SaveCrop(const std::string& aDir);
        void                        SaveMask(const std::string& aDir);
        void                        SaveSample(const std::string& aDir);
        void                        SaveTmp(tIm& aTmp, const std::string& aDir, const std::string& aPref);

    private:
        //----- members
        cPixBox<2>              mExtent;//-> input image CdT extent formed by predicted corners (input image coordinates)
        tIm                     mIm;    //-> input image
        tDIm*                   mDIm;   //-> input image data
        tIm                     mCrop;  //-> croped from input image
        tIm                     mRef;   //-> MMVII generated CdT
        tIm                     mSamp;  //-> CdT sampled from CdT2Im mapping
        tIm                     mMask;  //-> bbox of CdT formed by predicted corners (local coordinates, MaskIn/MaskOut)
        tU_INT1                 mVisib; //-> visibility score
        cAff2D_r                mRef2Im;
        cAff2D_r                mIm2Ref;
        std::vector<cPt2dr>     mVBCenters;//-> white bit centers w.r.t Ref image coordinates
        std::vector<cPt2dr>     mVWCenters;//-> black bit centers w.r.t Ref image coordinates
        cRansacSol              mTFSm2Cr;

        //----- methods
        void            SaveIm(tDIm* aDIm, std::string aPath);
        std::string     Cont();//-> returns actual content of the object
        /*      M : mMask OK    E : mExtent OK      OK i.e : tIm.Sz() != cPt2di(1,1)
         *      C : mCrop OK    R : mRef    OK
         */
    };

    cPixBox<2>          BBox(std::vector<cPt2dr> aVPts, int aMin=0, int aMax=100000);
    std::vector<cPt2dr> Corners(const cPt2dr& aP0, const cPt2dr& aP1);//-> get corners of a rectangle formed by aP0/aP1

    cAff2D_r Descr2Aff(const cCdTDescr& aDes, cSensorCamPC* aCam);//-> convert description to 2d affinity
}
