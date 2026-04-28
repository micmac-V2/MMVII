#ifndef C_EPIPOLAR_RECTIFICATION_H
#define C_EPIPOLAR_RECTIFICATION_H

#include "cPolyXY_N.h"
#include "MMVII_AllClassDeclare.h"  // cPt2dr, cPt3dr, cPt2di, etc.

namespace MMVII {
// Forward declaration
class cSensorImage;

class cEpipolarModel
{
public:
    virtual cPt2dr ToEpipolar1(const cPt2dr& p) const = 0;
    virtual cPt2dr ToEpipolar2(const cPt2dr& p) const = 0;
    virtual cPt2dr FromEpipolar1(const cPt2dr& p) const = 0;
    virtual cPt2dr FromEpipolar2(const cPt2dr& p) const = 0;
    cPt2dr ToEpipolar(int k, const cPt2dr& p) const
    {
        return k == 1 ? ToEpipolar1(p) : ToEpipolar2(p);
    }
    cPt2dr FromEpipolar(int k, const cPt2dr& p) const
    {
        return k == 1 ? FromEpipolar1(p) : FromEpipolar2(p);
    }
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
    //  Rectification result
    // --------------------------------------------------------
    struct cEpipolarModel : public MMVII::cEpipolarModel
    {
        // --- Rotation parameters (one per image) ---
        cPt2dr mCenter1, mCenter2;   ///< centroids of the image point sets
        cPt2dr mDir1,    mDir2;      ///< unit epipolar direction per image

        // --- Forward polynomials Vk (image -> epipolar) ---
        cPolyXY_N_IdentityOnYAxis<double> mV1;
        cPolyXY_N<double>                 mV2;

        // --- Inverse polynomials Wk (epipolar -> image) ---
        cPolyXY_N<double> mW1;
        cPolyXY_N<double> mW2;

        cEpipolarModel(int aPolyDeg, int aPolyDegInv);

        /// F1(p) : image-1 pixel -> epipolar pixel
        cPt2dr ToEpipolar1   (const cPt2dr& p) const override;
        /// F2(p) : image-2 pixel -> epipolar pixel
        cPt2dr ToEpipolar2   (const cPt2dr& p) const override;

        /// G1(e) : epipolar pixel -> image-1 pixel
        cPt2dr FromEpipolar1(const cPt2dr& e) const override;
        /// G2(e) : epipolar pixel -> image-2 pixel
        cPt2dr FromEpipolar2(const cPt2dr& e) const override;

    private:
        /// (p - C) / D  (complex division = rotation)
        static cPt2dr ToRotatedFrame  (const cPt2dr& p,
                                       const cPt2dr& aCenter,
                                       const cPt2dr& aDir);
        /// q * D + C
        static cPt2dr FromRotatedFrame(const cPt2dr& q,
                                       const cPt2dr& aCenter,
                                       const cPt2dr& aDir);
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
    cEpipolarModel Compute() const;

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
            const std::vector<cEpiPair>&        aPairs,
            cPolyXY_N_IdentityOnYAxis<double>&  aV1,
            cPolyXY_N<double>&                  aV2) const;

    // ----------------------------------------------------------
    //  Estimate inverse polynomials W1, W2 (eq. 33-34).
    //
    //  System :  Wk( qk.x ,  Vk(qk) ) = qk.y
    // ----------------------------------------------------------
    void EstimateInversePolynomials(
            const std::vector<cEpiPair>& aPairs,
            const cPolyXY_N<double>&     aV1,
            const cPolyXY_N<double>&     aV2,
            cPolyXY_N<double>&           aW1,
            cPolyXY_N<double>&           aW2) const;

    void EstimateOneInversePolynomial(
            const std::vector<cEpiPair>& aPairs,
            const cPolyXY_N<double>&     aVk,
            cPolyXY_N<double>&           aWk,
            bool                         aUsePt1) const;

    // --------------------------------------------------------
    //  Members
    // --------------------------------------------------------
    const cSensorImage& mCam1;
    const cSensorImage& mCam2;
    cParams             mParams;
};


// ---------------------------------------------------------------------------
// Layout descriptor for a single rectified image
// ---------------------------------------------------------------------------

struct EpipolarFrame {
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

// ---------------------------------------------------------------------------
// Full resampling result
// ---------------------------------------------------------------------------

struct EpipolarImages {
    cDataGenUnTypedIm<2>* im1_rect;   ///< rectified image 1
    cDataGenUnTypedIm<2>* im2_rect;   ///< rectified image 2
    EpipolarFrame frame1;     ///< layout descriptor for image 1
    EpipolarFrame frame2;     ///< layout descriptor for image 2
};

// ---------------------------------------------------------------------------
// Framing computation  (steps 1-3)
// ---------------------------------------------------------------------------

/**
 * Computes the rectified footprint of one image by projecting its original
 * boundary into epipolar space.
 *
 * Each of the four edges is densely sampled (nSample points per side) to
 * capture the curvature of epipolar lines, which can be significant for
 * pushbroom images.
 *
 * @param m        epipolar model
 * @param k        camera index (1 or 2)
 * @param ImSz    size of the original image (pixels)
 * @param nSample  number of sample points per edge (default 200)
 */
EpipolarFrame computeFrame(const cEpipolarModel& m, int k,
                           const cPt2di& ImSz,
                           int nSample = 200);

/**
 * Computes the common Y framing and the final dimensions of both rectified
 * images.
 *
 * Guarantees that the same row index v in im1_rect and im2_rect corresponds
 * to the same epipolar line.
 *
 * @param m              epipolar model
 * @param Im1Sz          dimensions of image 1
 * @param Im2Sz          dimensions of image 2
 * @param[out] f1, f2    computed layout descriptor for each image
 * @param margin         pixel margin added around the footprint (default 2)
 * @param nSample        edge sampling density (default 200)
 */
void computeCommonFraming(const cEpipolarModel& m,
                          const cPt2di& Im1Sz,
                          const cPt2di& Im2Sz,
                          EpipolarFrame& f1, EpipolarFrame& f2,
                          int margin  = 2,
                          int nSample = 200);

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
cDataGenUnTypedIm<2>* resampleToEpipolar(const cEpipolarModel& m,
                                         int k,
                                         const cDataGenUnTypedIm<2>* Im,
                                         const EpipolarFrame& frame,
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

EpipolarImages generateEpipolarImages(
    const cEpipolarModel& m,
    const cDataGenUnTypedIm<2>* Im1,
    const cDataGenUnTypedIm<2>* Im2,
    const cInterpolator1D& anInterp,
    int   margin  = 2,
    int   nSample = 200,
    double defVal  = 0);

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
