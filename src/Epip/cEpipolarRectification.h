#ifndef C_EPIPOLAR_RECTIFICATION_H
#define C_EPIPOLAR_RECTIFICATION_H

#include "cPolyXY_N.h"
#include "MMVII_Mappings.h"
#include "MMVII_AllClassDeclare.h"  // cPt2dr, cPt3dr, cPt2di, etc.
#include "MMVII_Tpl_ElemStrToVal.h"
#include "MMVII_MeasuresIm.h"  // cSetHomogCpleIm
#include <optional>

namespace MMVII {

class cSensorImage;

// Common vertical interval for epipolar Image pair resamling
enum class eEpipFrm
{
    eIntersect,         // Contains only common parts
    eUnion,             // Contains all parts
    eImg_1,             // Frame height from Image 1
    eImg_2,             // Frame height from Image 2
    eNbVals
};



class cEpipolarMapping : public cDataInvertibleMapping<tREAL8,2>
{
public:
    cEpipolarMapping(const cPt2dr& aZInterval, tREAL8 aGridStep, int aNbStepX, int aNbStepY)
        : mZInterval(aZInterval), mGridStep(aGridStep), mNbStepX(aNbStepX), mNbStepY(aNbStepY) {}
    void SetEpipImFrame(const cRect2& aFrame) { mEpipImFrame = aFrame;}
    cRect2 EpipFrame() const { return mEpipImFrame; }
    cPt2di EpipImSz() const { return mEpipImFrame.Sz(); }
    cPt2dr ZInterval() const { return mZInterval; }
    /// GenerateData's XY grid : pixel step and step count per axis
    tREAL8 GridStep() const { return mGridStep; }
    int NbStepX() const { return mNbStepX; }
    int NbStepY() const { return mNbStepY; }
protected:
    cRect2 mEpipImFrame{cPt2di{0,0},cPt2di{0,0},true}; ///< frame in epipolar space (for resampling)
    cPt2dr mZInterval;
    tREAL8 mGridStep;
    int    mNbStepX;
    int    mNbStepY;
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
                     cPt2dr aDir,
                     cPt2dr aZInterval,
                     tREAL8 aGridStep, int aNbStepX, int aNbStepY)
        : cEpipolarMapping(aZInterval, aGridStep, aNbStepX, aNbStepY)
        , mV(aV)
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
        int      mNbZLevels     = 3;      ///< number of altitude sampling levels
        eEpipFrm mEpipFrm       = eEpipFrm::eIntersect; ///< Framing type for epipolar images (Resmampling)
        int      mMargin        = 2;      ///< Margin in pixels for epipolar image framing (Resampling)
        std::optional<cPt2dr> mZIntv = std::nullopt; ///< Override Z interval (Zmin,Zmax); mandatory if sensor has none
        std::optional<cSetHomogCpleIm> mHomolPts = std::nullopt; ///< Tie points to infer Zmin/Zmax from (alternative to mZIntv, lower priority)
        tREAL8   mTiePMaxRes     = 2.0;   ///< Max triangulation residual (px) to keep a tie point for Z inference
        tREAL8   mZMargin        = 0.10;  ///< Relative margin added around the Z envelope inferred from mHomolPts
        tREAL8   mTiePMinNbRatio = 0.04;  ///< Min kept-tie-point count = max(mTiePMinNbFloor, mTiePMinNbRatio*sqrt(W*H))
        int      mTiePMinNbFloor = 25;    ///< Absolute floor for the min kept-tie-point count above
        bool     mNoWarnings     = false; ///< Don't generate warnigs: used by Bench
        size_t   mMinNbPairs     = 50;    ///< Min accepted pairs required in each of the train/test pools, per master-camera direction
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
    /// Independent (held-out) residuals : mean square of V1(q1)-V2(q2), W1(...)-q1.y, W2(...)-q2.y on the test pool
    double V1V2VarIndep() const { return mV1V2VarIndep; }
    double W1VarIndep() const { return mW1VarIndep; }
    double W2VarIndep() const { return mW2VarIndep; }
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
    //  Generate H-compatible pairs (Algorithm 2 of the paper), split into a
    //  train pool (fits V1/V2/W1/W2) and a test pool (EstimateIndepResiduals).
    //  Outputs: aOutPairsTrain/Test pairs, aOutCenterM centroid, aOutDirS direction.
    // ----------------------------------------------------------
    void GenerateData(const cSensorImage &aCamM, const cSensorImage &aCamS,
                      std::vector<cEpiPair> &aOutPairsTrain,
                      std::vector<cEpiPair> &aOutPairsTest,
                      cPt2dr &aOutCenterM,
                      cPt2dr &aOutDirS, cPt2dr &aZInterval,
                      tREAL8 &aOutGridStep, int &aOutNbStepX, int &aOutNbStepY) const;

    // ----------------------------------------------------------
    //  Independent residuals of the fitted V1,V2,W1,W2 on the held-out test pairs.
    // ----------------------------------------------------------
    void EstimateIndepResiduals(
            const std::vector<cEpiPair>& aPairsTest,
            const cPolyXY_Nd& aV1, const cPolyXY_Nd& aV2,
            const cPolyXY_Nd& aW1, const cPolyXY_Nd& aW2);

    // ----------------------------------------------------------
    //  Resolve the Z interval for a master camera : mZIntv > tie-point-derived
    //  > aCamM's own native interval.
    // ----------------------------------------------------------
    cPt2dr EffectiveZInterval(const cSensorImage & aCamM) const;

    // ----------------------------------------------------------
    //  Z envelope of mParams.mHomolPts, triangulated and filtered by
    //  residual, plus mZMargin. Memoized.
    // ----------------------------------------------------------
    cPt2dr ZIntervalFromHomolPts() const;
    mutable std::optional<cPt2dr> mCachedHomolZIntv;

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
    int mNbPairs12 = 0; ///< number of H-compatible training pairs from I1 to I2 (for info only)
    int mNbPairs21 = 0; ///< number of H-compatible training pairs from I2 to I1 (for info only)
    double mV1V2Var = 0.0;
    double mW1Var = 0.0;
    double mW2Var = 0.0;
    double mV1V2VarIndep = 0.0;
    double mW1VarIndep = 0.0;
    double mW2VarIndep = 0.0;
};


} // namespace MMVII

#endif // C_EPIPOLAR_RECTIFICATION_H
