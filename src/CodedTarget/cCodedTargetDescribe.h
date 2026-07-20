#ifndef CCODEDTARGETDESCRIBE_H
#define CCODEDTARGETDESCRIBE_H

#endif // CCODEDTARGETDESCRIBE_H

//#include "MMVII_Sensor.h"
#include "CodedTarget.h"
#include "MMVII_PCSens.h"
#include "MMVII_Geom3D.h"

namespace MMVII
{

extern const std::vector<cPt2dr> eSqCorners;//-> square corners coordinates from square centre
extern const tU_INT1 MaskOutV, MaskInV;//-> Val(aPix) = MaskOutV i.e aPix is out of the mask area

struct cCdTDetec;
struct cRansacSol;
struct cLS10PSol;
struct cCBParams;
class cAugCdt;
class cCdTDiscr;

typedef cSegment<tREAL8,3> tSeg3dr;
typedef cIm2D<tU_INT1>      tIm;
typedef cDataIm2D<tU_INT1>  tDIm;

/*!
     * @brief RansacATF : computes Affine Transfert Function between 2 images
     * @param aVBPts    : subset of classified points as black values
     * @param aVWPts    : subset of classified points as white values
     * @param aDIm1     : input image data
     * @param aDIm2     : output image data
     * @param aDMask    : mask image (no residual computation on points in the mask)
     * @param aGVT      : grey value treshold for i/o classification
     * @param aMaskOutV : value of out of the mask pixels
     * @param aIt       : number of iterations
     * @param aRDist    : minimum grey level distance between two w/b subset points
     * @return
     */

cRansacSol RansacATF(std::vector<cPt2dr> aVBPts, std::vector<cPt2dr> aVWPts, tDIm* aDIm1, tDIm* aDIm2,
                     tDIm* aDMask=nullptr, tU_INT1 aGVT=20, tU_INT1 aMaskOutV=MaskOutV, int aIt=200, int aRDist=50);

/*!
     * @brief The cRansacSol class stores RansacATF solution
     */

struct cRansacSol
{
    cRansacSol();
    cRansacSol(cPt2dr aSol);
    cRansacSol(cPt2dr aSol, cStdStatRes aL1Score, tIm aInlMask);
    cPt2dr      mSol;       //-> best fit solution
    cStdStatRes mL1Score;   //-> vector of inliers residuals
    tIm         mInlMask;   //-> inliers binary mask
};

struct cCBParams
{
    cCBParams();
    cCBParams(cPt2dr aC, tREAL8 aMLA, tREAL8 aLmA, tREAL8 aBCD);
    cPt2dr mC;//-> center
    tREAL8 mLMA;//-> lenght of semi-Major axis
    tREAL8 mLmA;//-> lenght of semi-minor axis
    tREAL8 mBCD;//-> bit center distance
};

/*!
     * @brief The cLS10PSol class stores cLS10PSys solution
     */

struct cLS10PSol
{
    cLS10PSol();
    cLS10PSol(cDenseVect<tREAL8> aVParams);
    std::vector<tREAL8> mVP;
    tAff2Dr            mAff;
};

/*!
     * @brief The cLS10PSys class represents the 10 parameters LS system to
     * and solve when computing Samp2Cr mapping
     */

class cLS10PSys
{
public:
    cLS10PSys(tDIm* aDIm1, tDIm* aDIm2, tDIm* aDInlMask, tU_INT1 aMaskInV=MaskInV);

    //-----methods
    void        Build();
    cLS10PSol   Solve();

private:
    //----- members
    tDIm*               mDIm1;
    tDIm*               mDIm2;
    tDIm*               mDInlMask;
    tU_INT1             mMaskInV;
    cLeasSqtAA<tREAL8>  mSys;
};

/*!
 * @brief The cExtract class stores CdT extractions
 */
struct cExtract
{
    cExtract(const cSensorCamPC* aCam, cSaveExtrEllipe aEll);
    const cSensorCamPC* mCam;
    const cSaveExtrEllipe mEll;
};

/*!
 * @brief The cAugCdt class performs CdT augmentation
 */
class cAugCdt
{
public:
    cAugCdt(std::string aName, std::shared_ptr<cFullSpecifTarget> aFSpec);
    cAugCdt();
    void                Spatialize(tREAL8 aGndInterTol=1e-2);
    cCdTDiscr           Discretize(cSensorCamPC* aCam, bool &isIn) const;
    void                AddExtract(cExtract aExt);
    tU_INT1             NbExtracts() const;
    void                AddData(const cAuxAr2007& anAux);
    std::vector<cPt3dr> GndCorners() const;
    std::string         Show() const;
    bool operator       <(const cAugCdt& aAug) const;
    std::string             mName;
    bool                    mOKAug;
    bool                    mOKInter;
    tREAL8                  m3DPrec;
    cSimilitud3D<tREAL8>    mRef2Gnd;

private:
    std::vector<cPt2dr>                 Corners() const;
    cPt2dr                              mCenter;
    std::shared_ptr<cFullSpecifTarget>  mFSpec;
    std::vector<cExtract>               mVExtracts;
};

void AddData(const cAuxAr2007& anAux, cAugCdt& anEx);

class cSetOfAugCdt
{
public:
    cSetOfAugCdt();
    void Add(cAugCdt aCdt);
    bool HasCdtName(std::string);
    void AddData(const cAuxAr2007& anAux);
    cAugCdt* AugOfName(std::string aName);
    std::vector<cAugCdt>& Cdts();
    static std::string NameFile(const cPhotogrammetricProject& aPhProj, bool aInput);
private:
    std::vector<cAugCdt> mVCdt;
};

void AddData(const cAuxAr2007& anAux, cSetOfAugCdt& aSet);


/*!
* @brief The cCdTDiscr class store and handles all about coded target discretisation
* for a specific image
*/

class cCdTDiscr
{
public:

    cCdTDiscr(const std::string& aName, const std::string& aImName, tAff2Dr aAff, bool aFromAug=false);
    cCdTDiscr();
    //----- members
    const std::string             mName;//-> CdT name
    const std::string             mImName;//-> name of the image in which CdT is discretized
    const bool                    mFromAug;//-> has it been created using network augmentation

    //----- methods
    bool                        ReqCt(std::string aRCt);//-> request content as a prerequisites check

    void                        Sample();//-> creates mSamp as a sampling of mCdT within mExtent
    void                        RansacTFSm2Cr();//-> computes mTFSm2Cr as a transfert function from mSamp to mCrop
    void                        VisuRansacTFSm2Cr(const std::string& aDir);//-> visu version of RansacTFSm2Cr : saves lots of visualisations
    void                        LS10ParamSm2Cr();//-> computes mSm2Cr mapping as a 10 parameters model
    void                        VisuLS10ParamSm2Cr(const std::string &aDir);//-> computes mSm2Cr mapping as a 10 parameters model
    tIm                         MaskInCB(bool ext=false);//-> computes Mask for pixels that are in checkboard pattern
    tIm                         MaskInCt(int aD);

    cPt2dr                      Ref2Im(cPt2dr aPt, bool inv=false);//-> get input image coordinate from a point of MMVII generated CdT image
    std::vector<cPt2dr>         VRef2Im(std::vector<cPt2dr> aVPts, bool inv=false);//-> Ref2Im for a point vector
    cPt2dr                      Ref2Samp(cPt2dr aP, bool inv=false);//-> convert aP w.r.t ref CdT coordinates to samp CdT coordinates

    void                        SetRef2Im(tAff2Dr aRef2Im);//-> set reference CdT image to input image mapping
    void                        SetExtent(cRect2 aExt);//-> set rectangular extent w.r.t input image coordinates
    void                        SetMask(tIm aMask);//-> set binary i/o CdT mask
    void                        SetInlMask(tIm aInlMask);//-> set binary i/o CdT inliers mask
    void                        SetRef(tIm aRef);//-> set reference image (MMVII generated)
    void                        SetCrop(tIm aCrop);//-> set cropped image
    void                        SetWBCenters(const std::vector<cPt2dr>& aV);//-> set w/b CdT bit centers from point vector
    void                        SetCB(std::unique_ptr<cFullSpecifTarget>& aFSpec);

    cRect2                      Extent();//-> rectangle extent of CdT in input image
    tIm&                        Mask();//-> i/o binary mask w.r.t samp image coordinates
    tIm&                        InlMask();//-> inliers mask for Samp to Crop mapping computation
    tIm&                        Crop();//-> input image cropped at mExtent bounds
    tIm&                        Ref();//-> referenced image (MMVII-generated)
    tIm&                        Samp();//-> sampled image
    cPt2dr                      SampC();//-> get sampled image center coordinate
    std::vector<cPt2dr>         SampBitC(char aCol);//-> returns 'w'-hite or 'b'-lack (aCol) bit centers w.r.t samp image coordinates
    cRansacSol                  TFSm2Cr();//-> get radiometric transfert function from sampled image to cropped image
    tIm                         ResTFSm2Cr();//-> computes residual image from sampled to cropped transfert function
    cLS10PSol                   LSSm2Cr();//-> get Samp to Crop mapping

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
    tIm                     mInlMask; //-> mask for inliers used for mapping LS adjustment
    tAff2Dr                 mRef2Glob;//-> reference (MMVII generated) CdT image to input image mapping
    std::vector<cPt2dr>     mVBCenters;//-> white bit centers w.r.t Ref image coordinates
    std::vector<cPt2dr>     mVWCenters;//-> black bit centers w.r.t Ref image coordinates
    cRansacSol              mTFSm2Cr;//-> radiometric transfert function from sampled image to cropped image computed by RansacATF
    cLS10PSol               mSm2Cr;//-> solution of LS 10 parameters solving
    cCBParams               mRefCB;//-> checkboard pattern parameters w.r.t Ref image coordinates

    //----- methods
    void            SaveIm(tDIm* aDIm, std::string aDir, std::string aPref);//-> generic image saving
    std::string     Cont();//-> returns actual content of the object in order to know what you can do
    /*      M : mMask OK    E : mExtent OK      OK i.e : tIm.Sz() != cPt2di(1,1)
         *      C : mCrop OK    R : mRef    OK
         */
    bool            InsideCB(cPt2dr aP, tREAL8 ext=0);//-> checks if a point with ref CdT coordinates is inside checkboard
};

class cSetOfCdTDiscr
{
public:
    cSetOfCdTDiscr(std::string aImName);
    void Add(cCdTDiscr aDis);
    bool HasCdTName(std::string);
    std::vector<cCdTDiscr> CdTDiscretizations() const;

private:
    std::vector<cCdTDiscr> mVDiscr;
    std::string mImName;
};

}
