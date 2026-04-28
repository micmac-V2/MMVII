#ifndef C_POLY_XY_N_H
#define C_POLY_XY_N_H

#include <vector>
#include <cassert>
#include <cmath>
#include <MMVII_Matrix.h>

// ============================================================
//  2D polynomial in variables (X,Y) of degree <= Degree :
//
//  P(x,y) = sum_{a+b <= Degree}  C[a,b] * x^a * y^b
//
//  Monomial ordering (lexicographic on (a,b)) :
//    (0,0),(0,1),...,(0,d),
//    (1,0),(1,1),...,(1,d-1),
//    ...
//    (d,0)
// ============================================================

namespace MMVII {


template <typename T>
class cPolyXY_N
{
public:
    // --------------------------------------------------------
    //  Construction
    // --------------------------------------------------------

    explicit cPolyXY_N(int aDegree)
        : mDegree(aDegree)
        , mCoeffs(NbCoeffsForDegree(aDegree), static_cast<T>(0))
    {}

    // --------------------------------------------------------
    //  Static helpers
    // --------------------------------------------------------

    /// Total number of monomials for a 2-variable polynomial of degree d
    static int NbCoeffsForDegree(int d);

    // --------------------------------------------------------
    //  Accessors
    // --------------------------------------------------------

    int Degree()   const { return mDegree; }
    int NbCoeffs() const { return static_cast<int>(mCoeffs.size()); }

    /// Linear index of monomial x^a * y^b  (requires a+b <= Degree)
    int Index(int a, int b) const;

    T &Coeff(int a, int b);

    T Coeff(int a, int b) const;

    T &CoeffByIdx(int k);

    T CoeffByIdx(int k) const;

    // --------------------------------------------------------
    //  Evaluation  P(x, y)
    // --------------------------------------------------------

    T Eval(T x, T y) const;

    /// Evaluation from a 2D point  P(p.x(), p.y())
    template <typename Pt2d>
    T Eval(const Pt2d &p) const;

    // --------------------------------------------------------
    //  Basis (monomial) vector at (x,y) for use in LSQ systems.
    //  aBase(Index(a,b)) = x^a * y^b
    // --------------------------------------------------------

    cDenseVect<T> BasisVector(T x, T y) const;

    template <typename Pt2d>
    cDenseVect<T> BasisVector(const Pt2d& p) const;

    // --------------------------------------------------------
    //  Fill coefficients from a LSQ solution vector
    // --------------------------------------------------------

    void SetFromSolution(const cDenseVect<T> &aSol, int aOffset = 0);

protected:
    int            mDegree;
    std::vector<T> mCoeffs;
};

typedef cPolyXY_N<double> cPolyXY_Nd;

template <typename T> inline int cPolyXY_N<T>::NbCoeffsForDegree(int d)
{
    return (d + 1) * (d + 2) / 2;
}


template <typename T> inline int cPolyXY_N<T>::Index(int a, int b) const
{
    assert(a >= 0 && b >= 0 && a + b <= mDegree);  // TODOCM: Use MMVII_ASSERT
    int idx = 0;
    for (int aa = 0; aa < a; ++aa)
        idx += (mDegree - aa + 1);
    idx += b;
    return idx;
    // TODOCM:
 //   MMVII_INTERNAL_ASSERT_medium(a>=0 && b>=0 && a+b<=mDegree,"Bad usage of cPolyXY_N");
 //   return b + ((mDegree + 1) * (mDegree + 2) - (mDegree - a + 1) * (mDegree - a + 2)) / 2;
}

template <typename T> inline T &cPolyXY_N<T>::Coeff(int a, int b)
{
    return mCoeffs[Index(a, b)];
}

template <typename T> inline T cPolyXY_N<T>::Coeff(int a, int b) const
{
    return mCoeffs[Index(a, b)];
}

template <typename T> inline T &cPolyXY_N<T>::CoeffByIdx(int k)
{
    assert(k >= 0 && k < NbCoeffs());  // TODOCM: Use MMVII_ASSERT
    return mCoeffs[k];
}

template <typename T> inline T cPolyXY_N<T>::CoeffByIdx(int k) const
{
    assert(k >= 0 && k < NbCoeffs());  // TODOCM: Use MMVII_ASSERT
    return mCoeffs[k];
}


template <typename T> inline T cPolyXY_N<T>::Eval(T x, T y) const
{
    T val = static_cast<T>(0);
    T px = static_cast<T>(1); // x^a
    for (int a = 0; a <= mDegree; ++a, px *= x) {
        T py = static_cast<T>(1); // y^b
        for (int b = 0; b <= mDegree - a; ++b, py *= y)
            val += mCoeffs[Index(a, b)] * px * py;
    }
    return val;
}

template <typename T>
template <typename Pt2d>
inline T cPolyXY_N<T>::Eval(const Pt2d &p) const
{
    return Eval(static_cast<T>(p.x()), static_cast<T>(p.y()));
}

template<typename T>
inline cDenseVect<T> cPolyXY_N<T>::BasisVector(T x, T y) const
{
    const int n = NbCoeffs();
    cDenseVect<T> v(n);
    for (int k = 0; k < n; ++k)         // TODOCM: supprimer ces 2 lignes
        v(k) = static_cast<T>(0);

    T px = static_cast<T>(1);
    for (int a = 0; a <= mDegree; ++a, px *= x)
    {
        T py = static_cast<T>(1);
        for (int b = 0; b <= mDegree - a; ++b, py *= y)
            v(Index(a, b)) = px * py;       // TODOCM : utiliser un idx++ a la place de Index(a,b)
    }
    return v;
}

template<typename T>
template<typename Pt2d>
inline cDenseVect<T> cPolyXY_N<T>::BasisVector(const Pt2d &p) const
{
    return BasisVector(static_cast<T>(p.x()),
                       static_cast<T>(p.y()));
}

template <typename T>
inline void cPolyXY_N<T>::SetFromSolution(const cDenseVect<T> &aSol,
                                          int aOffset)
{
    for (int k = 0; k < NbCoeffs(); ++k)
        mCoeffs[k] = aSol(aOffset + k);
}




// ============================================================
//  Specialised polynomial whose behaviour along the Y-axis
//  (i.e. when x=0) is locked to the identity :
//
//       P(0, y) = y   for all y
//
//  This is achieved by fixing the "x=0" coefficients:
//       C[0, 1] = 1
//       C[0, b] = 0   for b != 1
//
//  Geometrically this means: the polynomial maps the Y-axis
//  onto itself without distortion, which is required to make
//  the epipolar rectification unique (see Theorem 2, eq. 22).
// ============================================================

template <typename T>
class cPolyXY_N_IdentityOnYAxis : public cPolyXY_N<T>
{
public:
    using Base = cPolyXY_N<T>;

    explicit cPolyXY_N_IdentityOnYAxis(int aDegree);

    // --------------------------------------------------------
    //  Lock C[0,b] coefficients so that P(0,y) = y
    // --------------------------------------------------------

    void LockYAxisToIdentity();

    // --------------------------------------------------------
    //  Returns true if monomial (a,b) is free
    //  (i.e. not locked by the Y-axis identity constraint)
    // --------------------------------------------------------

    static bool IsFreeCoeff(int a, int /*b*/);

    /// Number of free (unlocked) coefficients
    int NbFreeCoeffs() const;

    // --------------------------------------------------------
    //  Basis vector restricted to free coefficients only.
    //  The locked part (a=0) is handled as a known RHS offset
    //  in the LSQ system.
    //
    //  Companion method : LockedContribution(x,y)  gives the
    //  value of the locked part at (x,y), i.e. the value of
    //  P(x,y) due solely to the fixed coefficients.
    //  => In the LSQ equation  V1(q1)=V2(q2),  the locked
    //     part must be moved to the RHS.
    // --------------------------------------------------------

    cDenseVect<T> FreeBasisVector(T x, T y) const;

    template <typename Pt2d> cDenseVect<T>
    FreeBasisVector(const Pt2d &p) const;

    /// Value contributed by the locked (fixed) coefficients at (x,y).
    /// Since C[0,1]=1 and all other C[0,b]=0 :  locked part = y
    T LockedContribution(T /*x*/, T y) const;

    template <typename Pt2d>
    T LockedContribution(const Pt2d &p) const;

    // --------------------------------------------------------
    //  Fill only the FREE coefficients from a solution vector.
    //  The locked coefficients are left untouched.
    // --------------------------------------------------------

    void SetFreeCoeffsFromSolution(const cDenseVect<T> &aSol, int aOffset = 0);
};


template <typename T>
inline cPolyXY_N_IdentityOnYAxis<T>::cPolyXY_N_IdentityOnYAxis(int aDegree)
    : Base(aDegree)
{
    LockYAxisToIdentity();
}

template <typename T>
inline void cPolyXY_N_IdentityOnYAxis<T>::LockYAxisToIdentity()
{
    for (int b = 0; b <= Base::mDegree; ++b)
        Base::Coeff(0, b) = (b == 1) ? static_cast<T>(1) : static_cast<T>(0);
}

template <typename T>
inline bool cPolyXY_N_IdentityOnYAxis<T>::IsFreeCoeff(int a, int /*b*/)
{
    return a >= 1;
}

template <typename T>
inline int cPolyXY_N_IdentityOnYAxis<T>::NbFreeCoeffs() const
{
    return Base::NbCoeffs() - (Base::mDegree + 1);
}

template <typename T>
inline cDenseVect<T> cPolyXY_N_IdentityOnYAxis<T>::FreeBasisVector(T x,
                                                                   T y) const
{
    const int d = Base::mDegree;
    const int nf = NbFreeCoeffs();
    cDenseVect<T> v(nf);
    for (int k = 0; k < nf; ++k)        // TODOCM: supprimer ces 2 lignes
        v(k) = static_cast<T>(0);

    T px = static_cast<T>(1);
    int idx = 0;
    for (int a = 0; a <= d; ++a, px *= x) {
        T py = static_cast<T>(1);
        for (int b = 0; b <= d - a; ++b, py *= y) {
            if (IsFreeCoeff(a, b))
                v(idx++) = px * py;
        }
    }
    return v;
}

template <typename T>
template <typename Pt2d>
inline cDenseVect<T>
cPolyXY_N_IdentityOnYAxis<T>::FreeBasisVector(const Pt2d &p) const
{
    return FreeBasisVector(static_cast<T>(p.x()), static_cast<T>(p.y()));
}

template <typename T>
inline T cPolyXY_N_IdentityOnYAxis<T>::LockedContribution(T /*x*/, T y) const
{
    return y; // = C[0,1] * y  = 1 * y
}

template <typename T>
template <typename Pt2d>
inline T cPolyXY_N_IdentityOnYAxis<T>::LockedContribution(const Pt2d &p) const
{
    return LockedContribution(static_cast<T>(p.x()), static_cast<T>(p.y()));
}

template <typename T>
inline void cPolyXY_N_IdentityOnYAxis<T>::SetFreeCoeffsFromSolution(
    const cDenseVect<T> &aSol, int aOffset)
{
    const int d = Base::mDegree;
    int idx = 0;
    for (int a = 0; a <= d; ++a)
        for (int b = 0; b <= d - a; ++b)
            if (IsFreeCoeff(a, b))
                Base::Coeff(a, b) = aSol(aOffset + idx++);
}

} // namespace MMVII
#endif // C_POLY_XY_N_H
