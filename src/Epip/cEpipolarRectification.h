#ifndef C_EPIPOLAR_RECTIFICATION_H
#define C_EPIPOLAR_RECTIFICATION_H

#include "cPolyXY_N.h"
#include "MMVII_Mappings.h"
#include "MMVII_AllClassDeclare.h"  // cPt2dr, cPt3dr, cPt2di, etc.

namespace MMVII {
// Forward declaration
class cSensorImage;

class cEpipolarMapping : public cDataInvertibleMapping<tREAL8,2>
{
public:
    cEpipolarMapping() {}
    cPt2dr ToEpipolar(const cPt2dr& aPt) const { return Value(aPt); }
    cPt2dr FromEpipolar(const cPt2dr& aPt) const { return Inverse(aPt); }
};


// ---------------------------------------------------------------------------
// Layout descriptor for a single rectified image
// ---------------------------------------------------------------------------

struct cEpipolarFrame {
    // Footprint in rectified space (coordinates phi_k applied to the original image)
    double xMin_rect = 0, xMax_rect = 0;   ///< X extent specific to image k
    double yMin_rect = 0, yMax_rect = 0;   ///< Y extent before intersection

    // Final common framing (Y shared between the two images)
    double yMin_common = 0, yMax_common = 0;

    // Integer offset applied to the output image.
    // Pixel (u, v) in the output <-> rectified coordinate (u + xOff, v + yOff)
    double xOff = 0;   ///< X offset, specific to image k
    double yOff = 0;   ///< Y offset, common to both images

    // Output image dimensions (pixels)
    int outSx = 0, outSy = 0;
};

class cEpipolarModel
{
public:
    virtual const cEpipolarMapping& EpipMap1() const = 0;
    virtual const cEpipolarMapping& EpipMap2() const = 0;
    cPt2dr ToEpipolar1(const cPt2dr& aPt) const { return EpipMap1().ToEpipolar(aPt); }
    cPt2dr FromEpipolar1(const cPt2dr& aPt) const { return EpipMap1().ToEpipolar(aPt); }
    cPt2dr ToEpipolar2(const cPt2dr& aPt) const { return EpipMap2().ToEpipolar(aPt); }
    cPt2dr FromEpipolar2(const cPt2dr& aPt) const { return EpipMap2().ToEpipolar(aPt); }

    void ComputeCommonFraming(const cPt2di& Im1Sz,
                              const cPt2di& Im2Sz,
                              int margin  = 2,
                              int nSample = 200);

    cDataGenUnTypedIm<2>* Resample1(const cDataGenUnTypedIm<2>* Im,
                                    const cInterpolator1D& anInterp,
                                    double defVal = 0.0)
    {
        return ResampleN(1,Im,anInterp,defVal);
    }

    cDataGenUnTypedIm<2>* Resample2(const cDataGenUnTypedIm<2>* Im,
                                    const cInterpolator1D& anInterp,
                                    double defVal = 0.0)
    {
        return ResampleN(2,Im,anInterp,defVal);
    }

    cDataGenUnTypedIm<2>* ResampleN(int Num,
                                    const cDataGenUnTypedIm<2>* Im,
                                    const cInterpolator1D& anInterp,
                                    double defVal = 0.0);


    cEpipolarFrame Frame1;
    cEpipolarFrame Frame2;

private:
    cEpipolarFrame ComputeFrame(const cEpipolarMapping& aEpipMap,
                                const cPt2di& ImSz,
                                int nSample = 200);

};


template<typename T>
#if __cplusplus >= 202002L
    requires std::derived_from<T, cEpipolarModelBase>
#endif
class cEpipolarModelTpl : public cEpipolarModel
{
public:
    typedef std::unique_ptr<T> Ptr_T;
    cEpipolarModelTpl(Ptr_T aEpipMap1, Ptr_T aEpipMap2)
        : aEpipMap1(std::move(aEpipMap1)),aEpipMap2(std::move(aEpipMap2))
    {}
    const cEpipolarMapping& EpipMap1() const override { return *aEpipMap1; }
    const cEpipolarMapping& EpipMap2() const override { return *aEpipMap2; }

private:
    Ptr_T aEpipMap1;
    Ptr_T aEpipMap2;
};



class cEpipPolyMapping: public cEpipolarMapping
{
public:
    cEpipPolyMapping(const cPolyXY_Nd& aV,
const cPolyXY_Nd& aW,
                     cPt2dr aCenter,
                     cPt2dr aDir)
        : mV(aV)
        , mW(aW)
        , mCenter{aCenter}
        , mDir{aDir}
    {}

    cPt2dr Value(const cPt2dr& aPt) const override;
    cPt2dr Inverse(const cPt2dr& aPt) const override;

private:
    /// (p - C) / D  (complex division = rotation)
    cPt2dr ToRotatedFrame(const cPt2dr& p) const;

    /// q * D + C
    cPt2dr FromRotatedFrame(const cPt2dr& q) const;

    // --- Forward polynomial Vk (image -> epipolar) ---
    cPolyXY_Nd mV;
    // --- Inverse polynomial Wk (epipolar -> image) ---
    cPolyXY_Nd mW;
    cPt2dr mCenter;   ///< centroids of the image point sets
    cPt2dr mDir;      ///< unit epipolar direction per image
};


class cEpipPolyModel : public cEpipolarModelTpl<cEpipPolyMapping>
{
public:
    using cEpipolarModelTpl<cEpipPolyMapping>::cEpipolarModelTpl;
};


// ============================================================
//  Epipolar rectification for a generic stereo camera pair.
//
//  Reference:
//    Pierrot-Deseilligny & Rupnik,
//    "Epipolar rectification of a generic camera", 2021.
//
//  The method computes two mappings
//
//    Fk : Ik -> Ek  (image -> epipolar space)
//
//  such that  F1(P1).y == F2(P2).y  whenever P1 and P2 are
//  H-compatible (i.e. they could be the projection of the same
//  3D point).
//
//  Each mapping has the form (eq. 23) :
//    Fk(p) = ( p_rot.x ,  Vk(p_rot) )
//    where  p_rot = (p - Ck) / Dk   (complex division = rotation)
//
//  The inverse is (eq. 33) :
//    Gk(e) = Dk * (e.x, Wk(e)) + Ck
// ============================================================

class cEpipolarRectification
{

public:
    // --------------------------------------------------------
    //  User parameters
    // --------------------------------------------------------
    struct cParams
    {
        int    mPolyDegree    = 3;      ///< degree of V polynomials
        int    mPolyDegreeInv = 7;      ///< degree of inverse W polynomials
        int    mNbXYSteps     = 100;    ///< number of image grid sampling steps (X & Y)
        int    mNbZLevels     = 3;      ///< number of altitude sampling levels
        double mEpsMarginRel  = 0.05;   ///< Relative Margin for X,Y and Z sampling
    };

    // --------------------------------------------------------
    //  Constructor
    // --------------------------------------------------------
    cEpipolarRectification(const cSensorImage& aCam1,
                           const cSensorImage& aCam2,
                           const cParams&      aParams);

    // --------------------------------------------------------
    //  Main entry point
    // --------------------------------------------------------
    cEpipPolyModel Compute() const;

private:
    // --------------------------------------------------------
    //  Private helper : one H-compatible pair in rotated coords
    // --------------------------------------------------------
    struct cEpiPair
    {
        cPt2dr mP1;   ///< rotated point in I1
        cPt2dr mP2;   ///< rotated point in I2
    };

    // ----------------------------------------------------------
    //  Generate H-compatible pairs (Algorithm 2 of the paper).
    //
    //  aCamM = master camera, aCamS = slave camera.
    //
    //  Outputs:
    //    aOutPairs : list of (masterPt, slavePt) pairs
    //    aOutCenterM : centroid of master image points
    //    aOutDirS    : average epipolar direction in slave image
    // ----------------------------------------------------------
    void GenerateData(const cSensorImage &aCamM, const cSensorImage &aCamS,
                      std::vector<cEpiPair> &aOutPairs, cPt2dr &aOutCenterM,
                      cPt2dr &aOutDirS) const;

    // ----------------------------------------------------------
    //  Estimate forward polynomials V1 (with Y-axis identity
    //  constraint) and V2 (unconstrained).
    //
    //  System (eq. 24) :  V1(q1) = V2(q2)
    // ----------------------------------------------------------
    void EstimateForwardPolynomials(
            const std::vector<cEpiPair>& aPairs,
            cPolyXY_Nd&           aV1,
            cPolyXY_Nd&           aV2) const;

    // ----------------------------------------------------------
    //  Estimate inverse polynomials W1, W2 (eq. 33-34).
    //
    //  System :  Wk( qk.x ,  Vk(qk) ) = qk.y
    // ----------------------------------------------------------
    enum class UseFromPair{PT1,PT2};
    void EstimateInversePolynomial(
            const std::vector<cEpiPair>& aPairs,
            const cPolyXY_Nd&     aVk,
            cPolyXY_Nd&           aWk,
            UseFromPair                  aUsePt) const;

    // --------------------------------------------------------
    //  Members
    // --------------------------------------------------------
    const cSensorImage& mCam1;
    const cSensorImage& mCam2;
    cParams             mParams;
};


// ---------------------------------------------------------------------------
// Resampling of one image into epipolar geometry  (step 4)
// ---------------------------------------------------------------------------

/**
 * Fills a rectified image from the original image.
 *
 * For each pixel (u, v) of the output image:
 *   rectified coordinate  e = (frame.xOff + u,  frame.yOff + v)
 *   original point        p = phi_k^{-1}(e)
 *   value                 outImg(u,v) = inImg.GetVBL(p)
 *
 * @param m       epipolar model
 * @param k       image index (1 or 2)
 * @param inImg   original image (accessed via DIm().GetVBL)
 * @param frame   layout computed by computeCommonFraming
 * @param defVal  default value for out-of-bounds pixels (default 0)
 */
cDataGenUnTypedIm<2>* resampleToEpipolar(
                                         const cEpipolarMapping &aEpipMap,
                                         const cDataGenUnTypedIm<2>* Im,
                                         const cEpipolarFrame& frame,
                                         const cInterpolator1D& anInterp,
                                         double defVal);



// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

/**
 * Generates both rectified images in the common epipolar geometry.
 *
 * Usage:
 *   EpipolarModel model = rectifier.compute();
 *   EpipolarImages result = generateEpipolarImages(model, im1, im2,
 *                                                  im1.DIm().Sx(), im1.DIm().Sy(),
 *                                                  im2.DIm().Sx(), im2.DIm().Sy());
 *   // result.im1_rect and result.im2_rect are the rectified images.
 *   // A pixel (u, v) in im1_rect is conjugate to (u', v) in im2_rect
 *   // for any u' (correspondence search is 1D along u').
 *
 * @param m               epipolar model (from EpipolarRectifier::compute())
 * @param Im1, Im2        original images
 * @param margin          pixel margin added around the footprint (default 2)
 * @param nSample         edge sampling density (default 200)
 * @param defVal          out-of-bounds fill value (default 0)
 */


#if 0
// ---------------------------------------------------------------------------
// Resampling quality check
// ---------------------------------------------------------------------------

/**
 * For a set of known conjugate point pairs (p1, p2), computes the residual
 * y-parallax in the rectified images.
 *
 * The V coordinate in rectified image k is:
 *   v_k = phi_k(p_k).y - frame_k.yOff
 *
 * The residual y-parallax is |v1 - v2|.
 *
 * @return max, mean, and RMS parallax statistics
 */
struct ResidualStats {
    double maxPar = 0, meanPar = 0, rmsPar = 0;
    size_t n = 0;
};

ResidualStats checkResiduals(const cEpipolarModel& m,
                             const EpipolarFrame& f1,
                             const EpipolarFrame& f2,
                             const std::vector<Pt2d>& pts1,
                             const std::vector<Pt2d>& pts2)
{
    if (pts1.size() != pts2.size())
        throw std::invalid_argument("checkResiduals: point vectors have different sizes");

    ResidualStats s;
    s.n = pts1.size();
    double sum = 0, sum2 = 0;

    for (size_t i = 0; i < s.n; ++i) {
        double v1  = toEpipolar(m, 1, pts1[i]).y - f1.yOff;
        double v2  = toEpipolar(m, 2, pts2[i]).y - f2.yOff;
        double par = std::abs(v1 - v2);
        s.maxPar = std::max(s.maxPar, par);
        sum  += par;
        sum2 += par * par;
    }
    if (s.n > 0) {
        s.meanPar = sum  / s.n;
        s.rmsPar  = std::sqrt(sum2 / s.n);
    }
    return s;
}


#endif
} // namespace MMVII

#endif // C_EPIPOLAR_RECTIFICATION_H
