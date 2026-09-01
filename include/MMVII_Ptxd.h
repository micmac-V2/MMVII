#ifndef  _MMVII_Ptxd_H_
#define  _MMVII_Ptxd_H_

#include "MMVII_AllClassDeclare.h"
#include "MMVII_nums.h"
#include "MMVII_memory.h"

#include <array>

namespace MMVII
{


/** \file MMVII_Ptxd.h
    \brief Basic N-dimensionnal point facilities

   Don't know exactly where we go ...  Probably will contain :

      - point in N dim template
      - specialization to 1,2,3,4  dims
      - dynamic dyn
      - boxes on others "small" connected classes

*/


template <class Type,const int Dim> class cPtxd;

template <class Type> class cDenseMatrix;
template <class Type> class cDenseVect;


///  template class for Points


template <class Type,const int Dim> class cPtxd
{
    public :
       typedef typename  tNumTrait<Type>::tBig               tBigNum ;
       typedef cPtxd<Type,Dim>                               tPt;
       typedef Type                                          tEl;
       // To see later (C Meynard ?) why this create compile pb
       // typedef typename  tElemNumTrait<Type>::tFloatAssoc    tReal;
       //typedef cPtxd<tReal,Dim>                              tPtR;


       static const int TheDim = Dim;
       /// Maybe some function will require generic access to data
       Type * PtRawData() {return mCoords.data();}
       const Type * PtRawData() const {return mCoords.data();}

       /// Unsafe acces to generic data
       Type & UncheckAt(int aK)
       {
           return mCoords[aK];
       }

       /// (const variant) Unsafe acces to generic data
       const Type & UncheckAt(int aK)  const
       {
           return mCoords[aK];
       }

       /// Safe acces to generic data
       Type & operator[] (int aK)
       {
          MMVII_INTERNAL_ASSERT_tiny((aK>=0) && (aK<Dim),"Bad point access");
          return mCoords[aK];
       }
       /// (const variant) Safe acces to generic data
       const Type & operator[] (int aK)  const
       {
          MMVII_INTERNAL_ASSERT_tiny((aK>=0) && (aK<Dim),"Bad point access");
          return mCoords[aK];
       }

       /// Access to data with a compilation constant value:  aPT::get<n>()
       template <int K>
       const Type &get() const {
           static_assert(K>=0 && K<Dim,"Invalid coordinate index in cPtxd.get<>()");
           return mCoords[K];
       }

       template <int K>
       Type &get() {
           static_assert(K>=0 && K<Dim,"Invalid coordinate index in cPtxd.get<>()");
           return mCoords[K];
       }

       /// Some function requires default constructor (serialization ?)
       cPtxd() : mCoords{} // fill array with 0
       { }

       /* I would prefer not inline but : cannot make work explicit instance of a
          a specific method,  and explicit of the whole class create problem
         with the static asser ;
       */
       // static cPtxd<Type,Dim> T_PCste(const Type & aVal) ;

       /// Initialisation with constants
       static cPtxd<Type,Dim>  PCste(const Type & aVal)
       {
           cPtxd<Type,Dim> aRes;
           tNumTrait<Type>::AssertValueOk(aVal);
           for (int aK=0 ; aK<Dim; aK++)
                aRes.mCoords[aK]= aVal;
           return aRes;
       }
       /// All coord 0 exect aCoord that values aVal
       static cPtxd<Type,Dim>  P1Coord(size_t aCoord,const Type & aVal) ;
       /// Initialisation with nan value (to detect error asap)
       static cPtxd<Type,Dim>  Dummy();
       /// Initialisation from name "i..."  "-j..."    valide are "ijkl"
       static cPtxd<Type,Dim>  PFromCanonicalName(const std::string & aName,size_t & anIndex);
      
       /// Initialisation from PInt
       static cPtxd<Type,Dim>  FromPtInt(const  cPtxd<int,Dim> & aVal) ;
       /// Initialisation from PInt
       static cPtxd<Type,Dim>  FromPtR(const  cPtxd<tREAL8,Dim> & aVal) ;
       /// Initialisation random
       static cPtxd<Type,Dim>  PRand();
       /// Initialisation random
       static cPtxd<Type,Dim>  PRandC();
       /// Initialisation random VUnit
       static cPtxd<Type,Dim>  PRandUnit();
       /// Pt random in sphere
       static cPtxd<Type,Dim>  PRandInSphere();

       /// Initialisation random VUnit not too close to P
       static cPtxd<Type,Dim>  PRandUnitDiff(const cPtxd<Type,Dim>&,const Type &aDist = 1e-2);
       /// Initialisation random VUnit not too close to P or -P
       static cPtxd<Type,Dim>  PRandUnitNonAligned(const cPtxd<Type,Dim>&,const Type &aDist = 1e-2);

        static cPtxd<Type,Dim> Col(const cDenseMatrix<Type>&,int aCol);  ///< Init with colum of matrix
        static cPtxd<Type,Dim> Line(int aLine,const cDenseMatrix<Type>&); ///< Init with line of matrix
        static cPtxd<Type,Dim> FromVect(const cDenseVect<Type>&); ///< Init with line of matrix
        static cPtxd<Type,Dim> FromStdVector(const std::vector<Type>&); ///< Init with line of matrix

        // Single-coordinate constructor: explicit
        explicit cPtxd(const Type& aCoord);

        // Multi-coordinate constructor: non-explicit
        template <typename... Args,
                 std::enable_if_t<(sizeof...(Args) > 1), int> = 0>      // forbid instanciation with only 1 arg, force use of explicit constructor
        cPtxd(const Args&... aCoords);

        explicit cPtxd(const Type * aV)
        {
           AssertTabValueOk(aV,size_t(Dim));
           MemCopy(&mCoords[0],aV,Dim);
        }

       // return Pt ot -Pt, the one "best" oriented with this
        tPt OrientInSameDir(const tPt& aPt) const;

        inline Type & x()             {return get<0>();}
        inline const Type & x() const {return get<0>();}

        inline Type & y()             {return get<1>();}
        inline const Type & y() const {return get<1>();}

        inline Type & z()             {return get<2>();}
        inline const Type & z() const {return get<2>();}

        inline Type & t()             {return get<3>();}
        inline const Type & t() const {return get<3>();}

        inline Type & u()             {return get<4>();}
        inline const Type & u() const {return get<4>();}

        cPtxd& operator+=(const cPtxd& aP);
        cPtxd& operator-=(const cPtxd& aP);
        cPtxd& operator*=(const Type& aVal);
        cPtxd& operator/=(const Type& aVal);

        cDenseVect<Type> ToVect() const; ///< conversion
        std::vector<Type> ToStdVector() const; ///< conversion
        void PushInStdVector(std::vector<Type> &) const;

        tBigNum  MinSqN2(const std::vector<tPt> &,bool SVP=false) const; ///< if SVP & empty return 0

        /// Used for "generik" object that must describes its box
        cTplBox<Type,Dim>  GetBoxEnglob() const;
        bool  InfEqDist(const tPt &,tREAL8) const;

        bool IsValid() const {
            for (int aK=0 ; aK<Dim; aK++)
                if  (!tNumTrait<Type>::ValueOk(mCoords[aK]))
                    return false;
            return true;
        }

        struct cDimIter {
            class iterator
            {
            public:
                iterator(int v) : mV(v) {}
                int operator*() const { return mV; }
                iterator& operator++() {++mV;return *this;}
                bool operator!=(const iterator& other) const {return mV != other.mV; }
            private:
                int mV;
            };
            iterator begin() const { return iterator(0); }
            iterator end() const { return iterator(Dim); }
        };

        static cDimIter DimIter() { return cDimIter{}; }

    protected :
        std::array<Type,Dim> mCoords;
};

template <class T,const int Dim>  class  cNV<cPtxd<T,Dim> >
{
    public :
        static  cPtxd<T,Dim>V0(){return  cPtxd<T,Dim>::PCste(0);}
};


// Explicit single coord constructor
template <class Type, const int Dim>
cPtxd<Type, Dim>::cPtxd(const Type& aCoord) :
    mCoords{aCoord}
{
    static_assert(Dim == 1, "Wrong number of coordinates in cPtdx initializer");
    tNumTrait<Type>::AssertValueOk(aCoord);
}

// Non-explicit mulit coords constructor
template <class Type, const int Dim>
template <typename... Args, std::enable_if_t<(sizeof...(Args) > 1), int>>
cPtxd<Type, Dim>::cPtxd(const Args&... aCoords) :
    mCoords{static_cast<Type>(aCoords)...}
{
    static_assert(sizeof...(Args) == Dim, "Wrong number of coordinates in cPtdx initializer");
    (tNumTrait<Type>::AssertValueOk(aCoords), ...);
}

// Arithmetic assignment operators
template <class Type, const int Dim>
cPtxd<Type, Dim>& cPtxd<Type, Dim>::operator+=(const cPtxd& aP)
{
    for (auto aK : DimIter()) this->UncheckAt(aK) += aP.UncheckAt(aK);
    return *this;
}

template <class Type, const int Dim>
cPtxd<Type, Dim>& cPtxd<Type, Dim>::operator-=(const cPtxd& aP)
{
    for (auto aK : DimIter()) this->UncheckAt(aK) -= aP.UncheckAt(aK);
    return *this;
}

// *= by a scalar
template <class Type, const int Dim>
cPtxd<Type, Dim>& cPtxd<Type, Dim>::operator*=(const Type& aVal)
{
    for (auto aK : DimIter()) this->UncheckAt(aK) *= aVal;
    return *this;
}

// /= by a scalar
template <class Type, const int Dim>
cPtxd<Type, Dim>& cPtxd<Type, Dim>::operator/=(const Type& aVal)
{
    MMVII_INTERNAL_ASSERT_NotNul(aVal);
    for (auto aK : DimIter()) this->UncheckAt(aK) /= aVal;
    return *this;
}

// Arithmetic operators
// +
template <class Type, const int Dim>
cPtxd<Type, Dim> operator+(cPtxd<Type, Dim> aP1, const cPtxd<Type, Dim>& aP2)
{
    cPtxd<Type, Dim> res;
    for (auto aK : aP1.DimIter()) res.UncheckAt(aK) = aP1.UncheckAt(aK) + aP2.UncheckAt(aK);
    return res;
}

// - binary
template <class Type, const int Dim>
cPtxd<Type, Dim> operator-(cPtxd<Type, Dim> aP1, const cPtxd<Type, Dim>& aP2)
{
    cPtxd<Type, Dim> res;
    for (auto aK : aP1.DimIter()) res.UncheckAt(aK) = aP1.UncheckAt(aK) - aP2.UncheckAt(aK);
    return res;
}

// - unary
template <class Type, const int Dim>
cPtxd<Type, Dim> operator-(const cPtxd<Type, Dim>& aP)
{
    cPtxd<Type, Dim> res;
    for (auto aK : aP.DimIter()) res.UncheckAt(aK) = -aP.UncheckAt(aK);
    return res;
}

// * scalar
template <class Type, const int Dim>
cPtxd<Type, Dim> operator*(cPtxd<Type, Dim> aP, const Type& aVal)
{
    cPtxd<Type, Dim> res;
    for (auto aK : aP.DimIter()) res.UncheckAt(aK) = aP.UncheckAt(aK) * aVal;
    return res;
}

// scalar *
template <class Type, const int Dim>
cPtxd<Type, Dim> operator*(const Type& aVal, cPtxd<Type, Dim> aP)
{
    return aP * aVal;
}

//  /  scalar
template <class Type, const int Dim>
cPtxd<Type, Dim> operator/(cPtxd<Type, Dim> aP, const Type& aVal)
{
    MMVII_INTERNAL_ASSERT_NotNul(aVal);
    cPtxd<Type, Dim> res;
    for (auto aK : aP.DimIter()) res.UncheckAt(aK) = aP.UncheckAt(aK) / aVal;
    return res;
}


// Scalar product (dot)
template <class Type, const int Dim>
Type Scal(const cPtxd<Type, Dim>& aP1, const cPtxd<Type, Dim>& aP2)
{
    auto res  = aP1.template get<0>() * aP2.template get<0>();  // cFormmula has no useable default constructor (Type = cFormula)
    for (int aK=1; aK < Dim; ++aK)
        res = res + aP1.UncheckAt(aK) * aP2.UncheckAt(aK);  // cFormula don't implement +=, can't write res+= aP1[Ak] + aP2[AK)
    return res;
}

// Scalar product (dot) with Big type
template <class Type, const int Dim>
typename tNumTrait<Type>::tBig ScalBig(const cPtxd<Type, Dim>& aP1, const cPtxd<Type, Dim>& aP2)
{
    auto res = static_cast<typename tNumTrait<Type>::tBig>(aP1.template get<0>()) * static_cast<typename tNumTrait<Type>::tBig>(aP2.template get<0>());
    for (int aK=1; aK < Dim; ++aK)
        res = res + static_cast<typename tNumTrait<Type>::tBig>(aP1.UncheckAt(aK)) * static_cast<typename tNumTrait<Type>::tBig>(aP2.UncheckAt(aK));
    return res;
}

// element-wise multiplication
template <class Type, const int Dim>
cPtxd<Type, Dim> MulCByC(const cPtxd<Type, Dim>& aP1, const cPtxd<Type, Dim>& aP2)
{
    cPtxd<Type, Dim> res;
    for (auto aK : aP1.DimIter()) res.UncheckAt(aK) = aP1.UncheckAt(aK) * aP2.UncheckAt(aK);
    return res;
}

// element-wise division ( integer division for integer type ! See RDivCByC)
template <class Type, const int Dim>
cPtxd<Type, Dim> DivCByC(const cPtxd<Type, Dim>& aP1, const cPtxd<Type, Dim>& aP2)
{
    cPtxd<Type, Dim> res;
    for (auto aK : aP1.DimIter()) res.UncheckAt(aK) = aP1.UncheckAt(aK) / aP2.UncheckAt(aK);
    return res;
}

// squared norm
template <class Type, const int Dim>
typename tNumTrait<Type>::tBig SqN2(const cPtxd<Type,Dim> & aP)
{
    return ScalBig(aP,aP);
}

// Generic comparison, col by col
template <class Type,int N,class Pred>
bool AllCmp(const cPtxd<Type,N>& aP1, const cPtxd<Type,N>& aP2, Pred pred)
{
    for (auto aK : aP1.DimIter())
        if (!pred(aP1.UncheckAt(aK),aP2.UncheckAt(aK)))
            return false;
    return true;
}

// ==
template <class Type,int N>
bool operator==(const cPtxd<Type,N>& aP1,const cPtxd<Type,N>& aP2)
{
    return AllCmp(aP1,aP2,std::equal_to<>{});
}

// !=
template <class Type,int N>
bool operator!=(const cPtxd<Type,N>& aP1,const cPtxd<Type,N>& aP2)
{
    return ! ( aP1 == aP2);
}

// element-wise >
template <class Type,int N>
bool SupStr(const cPtxd<Type,N>& aP1,const cPtxd<Type,N>& aP2)
{
    return AllCmp(aP1,aP2,std::greater<>{});
}

// element-wise >=
template <class Type,int N>
bool SupEq(const cPtxd<Type,N>& aP1,const cPtxd<Type,N>& aP2)
{
    return AllCmp(aP1,aP2,std::greater_equal<>{});
}

// element-wise <
template <class Type,int N>
bool InfStr(const cPtxd<Type,N>& aP1,const cPtxd<Type,N>& aP2)
{
    return AllCmp(aP1,aP2,std::less<>{});
}

// element-wise <=
template <class Type,int N>
bool InfEq(const cPtxd<Type,N>& aP1,const cPtxd<Type,N>& aP2)
{
    return AllCmp(aP1,aP2,std::less_equal<>{});
}

// Create the neighboord, ie pixel not nul, with coord in [-1,0,1]  having a  number of value  !=0  <= to aNbVois
// If dim =2 aNbVois->1 create the 4 neigh, NbVois-> 2 create the 8 neigh
// If Dim=3   1-> 6  2->    3->26 ( 3^3 -1)
template <const int Dim>  const std::vector<cPtxd<int,Dim>> & AllocNeighbourhood(int aNbVois);

//  classical  8-Neighbourhood
const std::vector<cPt2di> & Alloc8Neighbourhood();
//  classical  4-Neighbourhood
const std::vector<cPt2di> & Alloc4Neighbourhood();

//  Create a tab where K entrie represent vectors having NormInf equal to K
//  !! =>  Entry go from 0 to aDistMax included
//  !! =>   the size can be larger (but obviously not smaller) than dist required, as function remumber previous calls ....
template <const int Dim>  const std::vector<std::vector<cPtxd<int,Dim>>> & TabGrowNeigh(int aDistMax);

/// Return pixel between two radius, the order make them as sparse as possible (slow method in N^3) => To implement ???? No longer know what I wanted to do ???
//std::vector<cPt2di> SparsedVectOfRadius(const double & aR0,const double & aR1); // > R0 et <= R1
/// Implemented
std::vector<cPt2di> SortedVectOfRadius(const double & aR0,const double & aR1,bool IsSym=false); // > R0 et <= R1

/// IsSym  means that there is only one out of 2 between -P and P
std::vector<cPt2di> VectOfRadius(const double & aR0,const double & aR1,bool IsSym=false) ;



/** "Strange" function, because require DimIn=DimOut, but sometime we need to do this cast,
    probably consequence of bad design ... */
// Not needed ? (CM)
// template <class Type,const int DimOut,const int DimIn> void CastDim(cPtxd<Type,DimOut> &,const cPtxd<Type,DimIn>);

template <class Type, int Dim>
inline bool IsNull (const cPtxd<Type,Dim> & aP)
{
    for (auto aK : aP.DimIter()) {
        if (aP.UncheckAt(aK) != Type{})
            return false;
    }
    return true;
}

template <class Type, int Dim>
inline bool IsNotNull (const cPtxd<Type,Dim> & aP) { return ! IsNull(aP);}

#if (The_MMVII_DebugLevel>=The_MMVII_DebugLevel_InternalError_tiny )
template <class Type, int Dim> inline void AssertNonNul(const cPtxd<Type,Dim> &aP1)
{
   MMVII_INTERNAL_ASSERT_tiny(IsNotNull(aP1),"Unexpected null point");
}
#else
#define AssertNonNul(aP) {}
#endif


/// Some time int division is not what is wanted !!
template <class T,const int Dim> inline cPtxd<double,Dim> RDivCByC(const cPtxd<T,Dim> & aP1,const cPtxd<T,Dim> & aP2)
{
   return DivCByC(ToR(aP1),ToR(aP2));
}


///  Norms on points
template <class T,const int Dim> double NormK(const cPtxd<T,Dim> & aP,double anExp) ;
template <class T,const int Dim> T Norm1(const cPtxd<T,Dim> & aP);
template <class T,const int Dim> T NormInf(const cPtxd<T,Dim> & aP);

// template <class T,const int Dim> typename tNumTrait<T>::tFloatAssoc Norm2(const cPtxd<T,Dim> & aP);
template <class T,const int Dim> typename tNumTrait<T>::tFloatAssoc Norm2(const cPtxd<T,Dim> & aP)
{
   return sqrt(SqN2(aP));
}


///  x*y*z ...
template <class T,const int Dim> typename tNumTrait<T>::tBig MulCoord(const cPtxd<T,Dim> & aP);

template <class T,const int Dim> T Cos(const cPtxd<T,Dim> &,const cPtxd<T,Dim> &);
template <class T,const int Dim> T CosWDef(const cPtxd<T,Dim> &,const cPtxd<T,Dim> &,const T&);
template <class T,const int Dim> T AbsAngle(const cPtxd<T,Dim> &,const cPtxd<T,Dim> &);
///  Trunk cos in [-1,1] if necessary
template <class T,const int Dim> T AbsAngleTrnk(const cPtxd<T,Dim> &,const cPtxd<T,Dim> &);
/// angle of line, in [0,PI/2]  , opposite direction ar considered equal
template <class T,const int Dim> T AbsLineAngleTrnk(const cPtxd<T,Dim> &,const cPtxd<T,Dim> &);

//  +- equiv to AbsLineAngleTrnk for small values, more accurate and faster (?)
template <class T,const int Dim> T DistDirLine(const cPtxd<T,Dim> &,const cPtxd<T,Dim> &,const T & aDef=2.0);


template <class T,const int Dim> T MinAbsCoord(const cPtxd<T,Dim> & aP);


/// Sort vector by norm, typically dont need to compute square root
template <class Type,const int Dim> bool CmpN2(const cPtxd<Type,Dim> &aP1,const  cPtxd<Type,Dim> & aP2)
{
    return SqN2(aP1) < SqN2(aP2);
}


/// PtSupEq   : smallest point being SupEq to
template <class Type,int Dim>
cPtxd<Type,Dim> PtSupEq(const cPtxd<Type,Dim>& aP1, const cPtxd<Type,Dim>& aP2)
{
    cPtxd<Type,Dim> res;
    for (auto aK : aP1.DimIter()) res.UncheckAt(aK) = std::max(aP1.UncheckAt(aK),aP2.UncheckAt(aK));
    return res;
}

/// PtInfEq   : bigeest point being InfEq to
template <class Type,int Dim>
cPtxd<Type,Dim> PtInfEq(const cPtxd<Type,Dim>& aP1, const cPtxd<Type,Dim>& aP2)
{
    cPtxd<Type,Dim> res;
    for (auto aK : aP1.DimIter()) res.UncheckAt(aK) = std::min(aP1.UncheckAt(aK),aP2.UncheckAt(aK));
    return res;
}

template <class TypePt> void SetSupEq(TypePt & aP1,const TypePt & aP2) {aP1 = PtSupEq(aP1,aP2);}
template <class TypePt> void SetInfEq(TypePt & aP1,const TypePt & aP2) {aP1 = PtInfEq(aP1,aP2);}

template <class Type> inline cPtxd<Type,2> Transp  (const cPtxd<Type,2> & aP) {return  cPtxd<Type,2>(aP.y(),aP.x());}

/**  PtInfSTr : bigets point beg=ing InfStr (definition valide for integer types)
  Warn non symetric function;  strictness is relative to P2, not P1 ;
     For Floating point its inf as usuasl
*/
template <class Type> inline Type  VInfStr(const Type & aV1,const  Type & aV2) {return std::min(aV1,aV2);}
template <> inline int  VInfStr(const int & aV1,const  int & aV2) {return std::min(aV1,aV2-1);}
template <> inline tINT8  VInfStr(const tINT8 & aV1,const  tINT8 & aV2) {return std::min(aV1,aV2-1);}

template <class Type> inline cPtxd<Type,1> PtInfStr  (const cPtxd<Type,1> & aP1,const cPtxd<Type,1> & aP2)
{ return cPtxd<Type,1> (VInfStr(aP1.x(),aP2.x()));}

// { return cPtxd<Type,1> (std::min(aP1.x(),aP2.x()-1)); }
template <class Type> inline cPtxd<Type,2> PtInfStr  (const cPtxd<Type,2> & aP1,const cPtxd<Type,2> & aP2)
{ return cPtxd<Type,2> (VInfStr(aP1.x(),aP2.x()),VInfStr(aP1.y(),aP2.y()));}
// { return cPtxd<Type,2> (std::min(aP1.x(),aP2.x()-1),std::min(aP1.y(),aP2.y()-1)); }
template <class Type> inline cPtxd<Type,3> PtInfStr  (const cPtxd<Type,3> & aP1,const cPtxd<Type,3> & aP2)
{ return cPtxd<Type,3> (std::min(aP1.x(),aP2.x()-1),std::min(aP1.y(),aP2.y()-1),std::min(aP1.z(),aP2.z()-1)); }
template <class Type> inline cPtxd<Type,4> PtInfStr  (const cPtxd<Type,4> & aP1,const cPtxd<Type,4> & aP2)
{ return cPtxd<Type,4> (std::min(aP1.x(),aP2.x()-1),std::min(aP1.y(),aP2.y()-1),std::min(aP1.z(),aP2.z()-1),std::min(aP1.t(),aP2.t()-1)); }
template <class Type> inline cPtxd<Type,5> PtInfStr  (const cPtxd<Type,5> & aP1,const cPtxd<Type,5> & aP2)
{ return cPtxd<Type,5> (std::min(aP1.x(),aP2.x()-1),std::min(aP1.y(),aP2.y()-1),std::min(aP1.z(),aP2.z()-1),std::min(aP1.t(),aP2.t()-1),std::min(aP1.u(),aP2.u()-1)); }


template<class T,const int Dim> cPtxd<T,Dim>  VUnit(const cPtxd<T,Dim> & aP);
template <class TypePt> std::pair<TypePt,TypePt> OrthogonalizePair(const TypePt & aP1,const TypePt & aP2);


///  operator <<
template <class Type,const int Dim> std::ostream & operator << (std::ostream & OS,const cPtxd<Type,Dim> &aP);



template <class T,const int Dim> inline double RatioMax(const cPtxd<T,Dim> & aP1,const cPtxd<T,Dim> & aP2)
{
   return NormInf(RDivCByC(aP1,aP2));
}

/// Absolute value of surface of the parallelograp OP1 OP2, defined for dim=2&3
template <class Type,const int Dim>  Type AbsSurfParalogram(const cPtxd<Type,Dim>&aP1,const cPtxd<Type,Dim>&aP2);



// cPt2dr operator / (const cPt2dr &aP1,const cPt2dr & aP2) {return (aP1*conj(aP)}

// Most frequent conversion
template<int Dim>
inline cPtxd<int,Dim> ToI(const cPtxd<double,Dim> & aP) {
    cPtxd<int,Dim> aRes;
    for (auto aK : aP.DimIter())
        aRes.UncheckAt(aK) = round_ni(aP.UncheckAt(aK));
    return aRes;
}


template <class To, class From, int Dim>
cPtxd<To,Dim> ToType(const cPtxd<From,Dim>& aP)
{
    cPtxd<To,Dim> res;
    for (auto aK : aP.DimIter()) res.UncheckAt(aK) = static_cast<To>(aP.UncheckAt(aK));
    return res;
}


template <class Type, int Dim>
cPtxd<double,Dim> ToR(const cPtxd<Type,Dim>& aP)
{
    return ToType<tREAL8>(aP);
}

template <class Type, int Dim>
cPtxd<float,Dim> ToF(const cPtxd<Type,Dim>& aP)
{
    return ToType<tREAL4>(aP);
}

template <class T,const int Dim> cPtxd<tREAL8,Dim> Centroid(const std::vector<cPtxd<T,Dim> > & aVPts);
template <class T,const int Dim> cPtxd<tREAL8,Dim> Centroid(const std::vector<cPtxd<T,Dim>>&,const std::vector<T>&);
template <class T,const int Dim> cPtxd<T,Dim> Centroid(T aW0,const cPtxd<T,Dim> & aP0,const cPtxd<T,Dim> & aP1);
template <class T,const int Dim> cPtxd<T,Dim> Centroid(T aW0,const cPtxd<T,Dim> & aP0,T aW1,const cPtxd<T,Dim> & aP1);

template <class tContPts>  class cComputeCentroids
{
    public :
       typedef  typename tContPts::value_type tPts;
       typedef  typename tPts::tEl            tEl;

       static tPts   MedianCentroids(const tContPts &);
       static tPts   LinearWeigtedCentroids(const tContPts &,const tPts & aP0,double aSigma,double aExpN2=1.0,double aErRej=1e10);
       static tREAL8 SigmaDist(const tContPts &,const tPts & aP0,double aProp);

       static tPts  StdRobustCentroid(const tContPts &,double aProp,int aNbIter);

       /// Median of dist to  MedianCentroids
       static tREAL8   MedianSigma(const tContPts &);

       static tPts   LinearWeigtedCentroids(const tContPts &,const std::vector<tEl>* =nullptr);

};


template <class Type,int Dim,int aKth> bool  CmpCoord(const cPtxd<Type,Dim> & aP1,const cPtxd<Type,Dim> & aP2)
{
   static_assert((aKth>=0) && (aKth<Dim),"CmpCoord");
   return aP1[aKth] < aP2[aKth];
}


/**  CByC operator, apply an operator coordinate by coordinate, first version with one points */
template <class Type,const int Dim,class TypeFctr>
cPtxd<Type,Dim> CByC1P
                (
                   const cPtxd<Type,Dim>  & aP1,
                   const TypeFctr &         aFctr
                )
{
    cPtxd<Type,Dim> aRes;
    for (int aK=0 ; aK<Dim ; aK++)
        aRes[aK] = aFctr(aP1[aK]);
    return aRes;
}

/** Idem CByC when we need to change & force to int the result */
template <class Type,const int Dim,class TypeFctr>
cPtxd<int,Dim>  ICByC1P
                (
                   const cPtxd<Type,Dim>  & aP1,
                   const TypeFctr &         aFctr
                )
{
    cPtxd<int,Dim> aRes;
    for (int aK=0 ; aK<Dim ; aK++)
        aRes[aK] = aFctr(aP1[aK]);
    return aRes;
}
template <class Type,const int Dim>  cPtxd<int,Dim> Pt_round_down(const cPtxd<Type,Dim>&  aP);
template <class Type,const int Dim>  cPtxd<int,Dim> Pt_round_up(const cPtxd<Type,Dim>&  aP);
template <class Type,const int Dim>  cPtxd<int,Dim> Pt_round_ni(const cPtxd<Type,Dim>&  aP);
// template <class Type,const int Dim>  cPtxd<Type,Dim> PCste(const Type & aVal);


/**  CByC version with 2 points */
template <class Type,const int Dim,class TypeFctr>
cPtxd<Type,Dim> CByC2P
                (
                   const cPtxd<Type,Dim>  & aP1,
                   const cPtxd<Type,Dim>  & aP2,
                   const TypeFctr &         aFctr
                )
{
    cPtxd<Type,Dim> aRes;
    for (int aK=0 ; aK<Dim ; aK++)
        aRes[aK] = aFctr(aP1[aK],aP2[aK]);
    return aRes;
}

/// Number of pixel in square window
int NbPixVign(const int & aVign);
/// Number of pixel in a non square window
template <const int Dim> int NbPixVign(const cPtxd<int,Dim> & aVign);


/// Order coordinate so that it can define a box
template <class Type,const int Dim> void MakeBox(cPtxd<Type,Dim> & aP0,cPtxd<Type,Dim> & aP1)
{
    for (int aK=0 ; aK<Dim ; aK++)
        OrderMinMax(aP0[aK],aP1[aK]);
}


//  === TABULATION OF NEIGHBOORING

extern cPt2di  TAB4Corner[4] ; ///< {{1,1},{-1,1},{-1,-1},{1,-1}};



/**  Class for box, they are template as typically :
       - double will be used in geometric indexes QdTree or tiling
       - int will be used in bitmap manipulation
*/
template <class Type,const int Dim>  class cTplBox
{
    public :
        
        typedef Type                             tNum ;
        typedef typename  tNumTrait<Type>::tBig  tBigNum ;
        typedef cTplBox<Type,Dim>                tBox;
        typedef cPtxd<Type,Dim>                  tPt;
        typedef cPtxd<tBigNum,Dim>               tBigPt;
        static constexpr int                     NbCorners = 1<<Dim;
        typedef tPt   tCorner[NbCorners];

        cTplBox(const tPt & aP0,const tPt & aP1,bool AllowEmpty=false);
        cTplBox(const tPt & aSz,bool AllowEmpty=false); // Create a box with origin in 0,0,..
        static tBox Empty();
        static tBox FromVect(const tPt * aBegin,const tPt * aEnd,bool AllowEmpty=false);
        static tBox FromVect(const std::vector<tPt> & aVecPt,bool AllowEmpty=false);
        static tBox CenteredBoxCste(Type);
        static tBox BigBox();

        cTplBox<tREAL8,Dim> ToR() const;
        cTplBox<tINT4,Dim>  ToI() const;
        

        void AddData(const  cAuxAr2007 & anAux);

        // tPt & P0() {return mP0;} !!!!! BUUUG : CANNOT LET MODIFY BECAUSE OTHER FIELD NOT UP TO DATE
        // tPt & P1() {return mP1;} !!!!! BUUUG : CANNOT LET MODIFY BECAUSE OTHER FIELD NOT UP TO DATE
        const tPt & P0() const {return mP0;} ///< Origin of object
        const tPt & P1() const {return mP1;} ///< End of object
        const tPt & Sz() const {return mSz;} ///< Size of object
        tPt  Middle () const; ///< Middle of box

         // SEE BELLOW, IF USE RESYNCRONIZE OBJECT AFTER  , GIVE COMPLICATED NAME ON PURPOSE
        tPt & P0ByRef() {return mP0;} ///< Origin of object
        tPt & P1ByRef() {return mP1;} ///< End of object

        const tBigNum & NbElem() const {return mNbElem;}  ///< Surface  / Volume

        //const tPt & SzCum() const; ///< Cumulated size, rather internal use

        // Boolean operators
           /// Specialistion 1D
        bool Inside(const tNum & aX) const
        {
           // static_assert(Dim==1,"Bas dim for integer access");
           return (aX>=mP0.x()) && (aX<mP1.x());
        }
        /// Return closest point inside the box
        tPt  Proj(const tPt & aP) const {return PtInfStr(PtSupEq(aP,mP0),mP1);}
        /// Are the two box equals
        bool operator == (const tBox & aR2) const ;
        /// Is  this included in aB
        bool  IncludedIn(const  tBox & aB)const;
        /// Sometime we need to represent the empty box explicitely
        bool IsEmpty() const;
        tBox   Translate(const tPt & aPt)const;

        tBox Sup(const tBox & aBox) const;
        tBox Inter(const tBox & aBox)const; ///< Intersction handle empty case
        tBox Dilate(const tPt & aPt)const;  ///< Dilatation, as in morpho math : mP0-P mP1+P
        tBox Dilate(const Type & aVal)const;  ///< Dilatation with constant coordinate
        tBox ScaleCentered(const Type & aVal)const;  ///< Dilatation with constant scaling

        Type Insideness(const tPt & aP) const;

        /// Assert that it is inside
        template <class TypeIndex> void AssertInside(const TypeIndex & aP) const
        {
             MMVII_INTERNAL_ASSERT_tiny(Inside(aP),"Point out of image");
        }
        void AssertSameArea(const tBox & aV) const; ///<  Assert object are identic
        void AssertSameSz(const   tBox & aV) const;   ///<  Check only size

           /// Is this point/pixel/voxel  inside
        bool Inside(const tPt & aP) const  {return SupEq(aP,mP0) && InfStr(aP,mP1);}
           /// Specialistion 1D

        //  ---  object generation inside box ----------------

        tPt  FromNormaliseCoord(const cPtxd<double,Dim> &) const;  ///< [0,1] * => Rect
        cPtxd<double,Dim> ToNormaliseCoord(const tPt & aP) const;  ///< Rect => [0,1] *

        static cPtxd<double,Dim>  RandomNormalised() ;     ///<  Random point in "hyper cube" [0,1] ^ Dim
        tPt   GeneratePointInside() const;   ///< Random point in integer rect
        tBox  GenerateRectInside(double aPowSize=1.0) const; ///< Hig Power generate "small" rect, never empty

        static void Corners(tCorner & aRes,const tPt &aP0,const tPt &aP1);
        void Corners(tCorner & aRes) const;


        Type DistMax2Corners(const tPt&) const;
        static size_t NbFlagCorner() ;
        static tPt  CornerOfFlag(size_t aFlag,const tPt &aP0,const tPt &aP1) ;
        tPt  CornerOfFlag(size_t aFlag) const;

    protected :
        tPt       mP0;         ///< "smallest"
        tPt       mP1;         ///< "highest"
        tPt       mSz;         ///<  Size
        tBigPt    mSzCum;      ///< Cumlated size : Cum[aK] = Cum[aK-1] * Sz[aK-1]
        tBigNum   mNbElem;     ///< Number of pixel = Cum[Dim-1]
    private :

};

/**  Assure that P0,P1 are non empty box, using a minimum changes */
template <const int Dim> void  MakeBoxNonEmptyWithMargin(cPtxd<tREAL8,Dim>&P0 , cPtxd<tREAL8,Dim> & aP1,
                                                         tREAL8 aStdMargin,tREAL8 aMarginSemiEmpty,tREAL8 aMarginEmpty);


// template <const int Dim>  cTplBox<tREAL8,Dim> ToR(const  cTplBox<int,Dim> & );
// template <const int Dim>  cTplBox<int,Dim> ToI(const  cTplBox<tREAL8,Dim> & );

/** Function computing corner of box, this one is specific to dim=1 because it respect
trigonometric order, a notion not generalisable */

template <class Type> void CornersTrigo(typename cTplBox<Type,2>::tCorner & aRes,const cTplBox<Type,2>&);

/*
typedef cTplBox<int,2>  cBox2di;
typedef cTplBox<double,2>  cBox2dr;
typedef cTplBox<int,3>  cBox3di;
typedef cTplBox<double,3>  cBox3dr;
*/



//cBox2dr ToR(const cBox2di &);  ///< Basic conversion
//cBox2di ToI(const cBox2dr &);  ///< Convert in englobing mode
cBox2dr operator * (const cBox2dr & aBox,double aScale); ///< just multiply each coord


// Is window inside the box
template <class Type> bool WindInside(const cBox2di & aBox,const cPtxd<Type,2> & aPt,const  cPt2di & aSzW);
template <class Type> bool WindInside(const cBox2di & aBox,const cPtxd<Type,2> & aPt,const  int & aSzW);
// Is window inside the box taking into account bilinear interpol ?
template <class Type> bool WindInside4BL(const cBox2di & aBox,const cPtxd<Type,2> & aPt,const  cPt2di & aSzW);


cBox2di DilateFromIntervPx(const cBox2di & aBox,int aDPx0,int aDPx1);


template <class Type,const int Dim> std::ostream & operator << (std::ostream & OS,const cTplBox<Type,Dim> &aBox)
{ return  OS << "{" << aBox.P0() <<   " :: " << aBox.P1()<< "}"; }

/**  Class for computing box of a set of points by iteratively adding them.
     Is ok with empty case (!= cTplBox)
     Can be converted to a "regular" box (cTplBox)
*/

template <class Type,const int Dim>  class cTplBoxOfPts
{
    public :
        typedef cPtxd<Type,Dim>                  tPt;

        cTplBoxOfPts();
        static cTplBoxOfPts FromVect(const tPt * aBegin,const tPt * aEnd);
        static cTplBoxOfPts FromVect(const std::vector<tPt> & aVecPt);

        int NbPts() const;  ///< Use to check acces that are forbidden when empty
        const tPt & P0() const;
        const tPt & P1() const;
        cTplBox<Type,Dim> CurBox(bool AllowEmpty=false) const;

        void Add(const tPt &);
        void AddData(const  cAuxAr2007 & anAux);
    private :
        int  mNbPts;  ///< Number of points, to check access
        tPt  mP0;
        tPt  mP1;
};

template <class Type,const int Dim> class cSegment
{
    public :
       typedef cPtxd<Type,Dim> tPt;
       cSegment(const tPt& aP1,const tPt& aP2);
       /// Estimate fonc linear, with gradient paral to tangent,  given value in P1 and P2, will be F(Q) =  R.first + R.second Q
       void CompileFoncLinear(Type & aVal,tPt & aVec,const Type  &aV1,const Type  & aV2) const;
       const tPt&  P1() const {return mVPts.at(0);}   ; ///< Accessor
       const tPt&  P2() const {return mVPts.at(1);}   ; ///< Accessor
       tPt&  P1()  {return mVPts.at(0);}   ; ///< Accessor
       tPt&  P2()  {return mVPts.at(1);}   ; ///< Accessor
       void  Swap(); ///< swap P1 P2 in place
       const std::vector<tPt> &  VPts() const;

       tPt  V12() const;   ///<  Vector  P1->P2
       tPt  Middle() const;  ///<  P middle
        /// Used for "generik" object that must describes its box
        cTplBox<Type,Dim>  GetBoxEnglob() const;
    protected :
       std::vector<tPt> mVPts;
       // tPt  mP1;
       // tPt  mP2;
};

template <class Type,const int Dim> class cSegmentCompiled : public cSegment<Type,Dim>
{
    public :
       typedef cPtxd<Type,Dim> tPt;
       cSegmentCompiled(const tPt& aP1,const tPt& aP2);
       cSegmentCompiled(const cSegment<Type,Dim>&);
       tPt  Proj(const tPt &) const;
       Type Dist(const tPt &) const; // dist to full line
       Type Abscissa(const tPt& aPt) const;
       tPt  PtOfAbscissa(const Type & anAbsc) const;

       const Type & N2 () const;
       const tPt  & Tgt() const;
    protected :
       Type    mN2;
       tPt     mTgt;
};

/// class for modelization of an affine space : a point + vectorials
template<const int Dim> class cAffineSpace
{
   public :
      typedef cPtxd<tREAL8,Dim> tPtR;
      typedef std::vector<tPtR> tVecSp;
      typedef cAffineSpace<Dim> tAffSp;

      const tPtR &   P0()    const {return mP0;}
      const tVecSp & VecSp() const {return mVecSp;}

      static tAffSp LstSqEstimate(const  std::vector<tPtR> &,int aSzSubs,const std::vector<tREAL8> * aVW=nullptr);

   private :
      tPtR    mP0;
      tVecSp  mVecSp;
};



};


template <class Type>
struct std::less<MMVII::cPtxd<Type,1>>
{
    std::size_t operator()(const MMVII::cPtxd<Type,1>& lhs, const MMVII::cPtxd<Type,1>& rhs) const
    {
        return lhs.x() < rhs.x();
    }
};

template <class Type>
struct std::less<MMVII::cPtxd<Type,2>>
{
    std::size_t operator()(const MMVII::cPtxd<Type,2>& lhs, const MMVII::cPtxd<Type,2>& rhs) const
    {
        return lhs.x() < rhs.x() ? true : lhs.y() < rhs.y();
    }
};

template <class Type>
struct std::less<MMVII::cPtxd<Type,3>>
{
    std::size_t operator()(const MMVII::cPtxd<Type,3>& lhs, const MMVII::cPtxd<Type,3>& rhs) const
    {
        return lhs.x() < rhs.x() ? true : (lhs.y() < rhs.y() ? true : lhs.z() < rhs.z());
    }
};

template <class Type>
struct std::less<MMVII::cPtxd<Type,4>>
{
    std::size_t operator()(const MMVII::cPtxd<Type,4>& lhs, const MMVII::cPtxd<Type,4>& rhs) const
    {
        return lhs.x() < rhs.x() ? true : (lhs.y() < rhs.y() ? true : (lhs.z() < rhs.z() ? true : lhs.t() < rhs.t()));
    }
};


#endif  //  _MMVII_Ptxd_H_
