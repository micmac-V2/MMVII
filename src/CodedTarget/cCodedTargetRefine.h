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

    tRect2          BBox(std::vector<cPt2dr> aVPts, int aMin=0, int aMax=100000);//-> computes bounding box from a point vector
    std::vector<cPt2dr> Corners(const cPt2dr& aP0, const cPt2dr& aP1);//-> gets corners of a rectangle formed by aP0/aP1


    class cPixBBox : public cPixBox<2>
    {
    public:
        typedef cPtxd<tREAL8,2> tPt;
        cPixBBox(const std::vector<tPt> aVPts);
        cPixBBox(cPixBox<2>);
    };

    template <class Type> class cOptCorrelThIm : public cDataMapping<tREAL8, 2, 1>
    {
    public:
        typedef  cDataIm2D<Type> tDIm;
        cOptCorrelThIm(tDIm& aTheorDIm, tDIm& aGlobDIm, cDataIm2D<tU_INT1>& aMaskDIm, cPixBox<2> aBBox);
        cPt1dr Value(const cPt2dr& aPt) const override;//-> correl score from aPt position
    private:
        tDIm& mThDIm;//-> theoretical image of the target generated from detected deformation
        tDIm& mGDIm;//-> global image
        tDIm& mDMask;//-> mask for correlation computation
        cPixBox<2> mBBox;//-> bbox of predicted target wrt global image
        cPt2dr mP0;//-> initial bbox up corner to compute translation
    };

    /**
     * @brief The cSampler class sample template image into a global image
     * according to a map and the position of the template in the global image
     */

    class cSampler
    {
    public:
        cSampler(tIm aInIm, cPixBox<2> aInBox, tAff2Dr aIn2Out);
        std::pair<tIm, cPt2di> Sample();
    private:
        tIm mIm;//-> original image to sample
        cPixBox<2> mBBox;//->
        tAff2Dr mMap;
    };

    /**
     * @brief The cHeuristSeeker class finds template image best position in a
     * global image using Heuristik optimization
     */

    class cHeuristSeeker
    {
    public:
        cHeuristSeeker(tIm& aCdtIm, tIm& aGlobIm);
        cPt2dr Seek(cPt2dr aP0);
    private:
        tIm& mTIm;
        tIm& mGIm;
    };

    /**
     * @brief The cPatch class represents a patch inside a global image
     */

    class cPatch
    {
    public:
        cPatch(tIm& aGlobIm, cPixBBox aGlobBBox);
        tIm Patch();
        void SaveIm(std::string aDir, bool isOk);
        cPixBox<2> BBox();
    private:
        tIm& mGIm;
        cPixBox<2> mBBox;
    };

    class cMaskO2I
    {
    public:
        cMaskO2I(const cPixBox<2>& aOBox, const cPixBox<2>& aIBox, const tAff2Dr& aO2IMap);
        void SaveAsIm(const std::string& aDir);
        tIm Im();
    private:
        const cPixBox<2>& mOBox;
        const cPixBox<2>& mIBox;
        const tAff2Dr& mO2IMap;
        tIm mIm;
    };

    class cBox2DTransfo
    {
    public:
        cBox2DTransfo(cPixBox<2>& aIBox, tAff2Dr aMap);
        cPixBox<2> Transfo();
    private:
        cPixBox<2>& mIBox;
        tAff2Dr mMap;
    };
}
