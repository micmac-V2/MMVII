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
     * @param aGVT      : grey value treshold for i/o classification
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
        tIm         mInlMask;//-> inliers binary mask
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
        const std::string             mName;//-> CdT name
        const std::string             mImName;//-> name of the image in which CdT is discretized

        //----- methods
        bool                        ReqCt(std::string aRCt);//-> request content as a prerequisites check

        void                        Sample();//-> creates mSamp as a sampling of mCdT within mExtent
        void                        RansacTFSm2Cr();//-> computes mTFSm2Cr as a transfert function from mSamp to mCrop
        void                        VisuRansacTFSm2Cr(const std::string& aDir);//-> visu version of RansacTFSm2Cr : saves lots of visualisations

        cPt2dr                      Ref2Im(cPt2dr aPt, bool inverse=false);//-> get input image coordinate from a point of MMVII generated CdT image
        std::vector<cPt2dr>         VRef2Im(std::vector<cPt2dr> aVPts, bool inverse=false);//-> Ref2Im for a point vector

        void                        SetRef2Im(cAff2D_r aRef2Im);//-> set reference CdT image to input image mapping
        void                        SetExtent(cRect2 aExt);//-> set rectangular extent w.r.t input image coordinates
        void                        SetMask(tIm& aMask);//-> set binary i/o CdT mask
        void                        SetRef(tIm& aRef);//-> set reference image (MMVII generated)
        void                        SetCrop(tIm& aCrop);//-> set cropped image
        void                        SetWBCenters(const std::vector<cPt2dr>& aV);//-> set w/b CdT bit centers from point vector

        cRect2                      Extent();//-> rectangle extent of CdT in input image
        tIm&                        Mask();//-> i/o binary mask w.r.t samp image coordinates
        tIm&                        Crop();//-> input image cropped at mExtent bounds
        tIm&                        Ref();//-> referenced image (MMVII-generated)
        tIm&                        Samp();//-> sampled image
        std::vector<cPt2dr>         SampC(char aC);//-> returns w/b bit centers w.r.t samp image coordinates
        cRansacSol                  TFSm2Cr();//-> get radiometric transfert function from sampled image to cropped image
        tIm                         ResTFSm2Cr();//-> computes residual image from sampled to cropped transfert function

        void                        SaveCrop(const std::string& aDir);
        void                        SaveMask(const std::string& aDir);
        void                        SaveSample(const std::string& aDir);
        void                        SaveRef(const std::string& aDir);
        void                        SaveTmp(tIm& aTmp, const std::string& aDir, const std::string& aPref);//-> save images which aren't class members

    private:
        //----- members
        cPixBox<2>              mExtent;//-> input image CdT extent formed by predicted corners (input image coordinates)
        tIm                     mIm;    //-> input image
        tDIm*                   mDIm;   //-> input image data
        tIm                     mCrop;  //-> croped from input image
        tIm                     mRef;   //-> MMVII generated CdT
        tIm                     mSamp;  //-> CdT sampled from CdT2Im mapping
        tIm                     mMask;  //-> bbox of CdT formed by predicted corners (local coordinates, MaskIn/MaskOut)
        cAff2D_r                mRef2Im;//-> reference (MMVII generated) CdT image to input image mapping
        std::vector<cPt2dr>     mVBCenters;//-> white bit centers w.r.t Ref image coordinates
        std::vector<cPt2dr>     mVWCenters;//-> black bit centers w.r.t Ref image coordinates
        cRansacSol              mTFSm2Cr;//-> radiometric transfert function from sampled image to cropped image computed by RansacATF

        //----- methods
        void            SaveIm(tDIm* aDIm, std::string aDir, std::string aPref);//-> generic image saving
        std::string     Cont();//-> returns actual content of the object in order to know what you can do
        /*      M : mMask OK    E : mExtent OK      OK i.e : tIm.Sz() != cPt2di(1,1)
         *      C : mCrop OK    R : mRef    OK
         */
    };

    cPixBox<2>          BBox(std::vector<cPt2dr> aVPts, int aMin=0, int aMax=100000);//-> compute bounding box from a point vector
    std::vector<cPt2dr> Corners(const cPt2dr& aP0, const cPt2dr& aP1);//-> get corners of a rectangle formed by aP0/aP1
    cAff2D_r Descr2Aff(const cCdTDescr& aDes, cSensorCamPC* aCam);//-> convert description to 2d affinity
}
