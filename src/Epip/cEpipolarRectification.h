#ifndef C_EPIPOLAR_RECTIFICATION_H
#define C_EPIPOLAR_RECTIFICATION_H

#include "cPolyXY_N.h"
#include "MMVII_Mappings.h"
#include "MMVII_AllClassDeclare.h"  // cPt2dr, cPt3dr, cPt2di, etc.

namespace MMVII {

class cSensorImage;

class cEpipolarMapping : public cDataInvertibleMapping<tREAL8,2>
{
public:
    cEpipolarMapping() {}
    void SetEpipImFrame(const cRect2& aFrame) { mEpipImFrame = aFrame;}
    cRect2 EpipFrame() const { return mEpipImFrame; }
    cPt2di EpipImSz() const { return mEpipImFrame.Sz(); }
protected:
    cRect2 mEpipImFrame{cPt2di{0,0},cPt2di{0,0},true}; ///< frame in epipolar space (for resampling)
};




class cEpipolarModel
{
public:
    virtual ~cEpipolarModel() {}
    const cEpipolarMapping& EpipMap1() const { return (const_cast<cEpipolarModel*>(this))->GetEpipMap1();}
    const cEpipolarMapping& EpipMap2() const { return (const_cast<cEpipolarModel*>(this))->GetEpipMap2();;}
    void ComputeCommonFraming(
        const cTplBox<tREAL8,2> aBox1,
        const cTplBox<tREAL8,2> aBox2,
        eEpipFrm aFrmType=eEpipFrm::eIntersect,
        int aMargin=0
    );
protected:
    virtual cEpipolarMapping& GetEpipMap1() = 0;
    virtual cEpipolarMapping& GetEpipMap2() = 0;

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
        : aEpipMap1(std::move(aEpipMap1)),aEpipMap2(std::move(aEpipMap2)  )
    {}

private:
    cEpipolarMapping& GetEpipMap1() override { return *aEpipMap1; }
    cEpipolarMapping& GetEpipMap2() override { return *aEpipMap2; }

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
        int      mPolyDegree    = 3;      ///< degree of V polynomials
        int      mPolyDegreeInv = 7;      ///< degree of inverse W polynomials
        int      mNbXYSteps     = 100;    ///< number of image grid sampling steps (X & Y)
        int      mNbZLevels     = 3;      ///< number of altitude sampling levels
        eEpipFrm mEpipFrm       = eEpipFrm::eIntersect; ///< Framing type for epipolar images (Resmampling)
        int      mMargin        = 2;      ///< Margin in pixels for epipolar image framing (Resampling)
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
    cEpipPolyModel Compute();

    int NbPairs12() const { return mNbPairs12; }
    int NbPairs21() const { return mNbPairs21; }
    double V1V2Var() const { return mV1V2Var; }
    double W1Var() const { return mW1Var; }
    double W2Var() const { return mW2Var; }
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
            cPolyXY_Nd&           aV2);

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
            UseFromPair                  aUsePt);

    // --------------------------------------------------------
    //  Members
    // --------------------------------------------------------
    const cSensorImage& mCam1;
    const cSensorImage& mCam2;
    cParams             mParams;
    int mNbPairs12 = 0; ///< number of H-compatible pairs from I1 to I2 (for info only)
    int mNbPairs21 = 0; ///< number of H-compatible pairs from I2 to I1 (for info only)
    double mV1V2Var = 0.0;
    double mW1Var = 0.0;
    double mW2Var = 0.0;
};


} // namespace MMVII

#endif // C_EPIPOLAR_RECTIFICATION_H
