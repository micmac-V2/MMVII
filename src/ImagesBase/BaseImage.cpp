#include "MMVII_Images.h"
#include "MMVII_Image2D.h"
#include "MMVII_Mappings.h"
#include <algorithm>
// #include <Eigen/Dense>

namespace MMVII
{

/* ========================== */
/*     cDataGenUnTypedIm      */
/* ========================== */


template <const int Dim> cDataGenUnTypedIm<Dim>::cDataGenUnTypedIm
                         (
                             const tPixI & aP0,
                             const tPixI & aP1
                         )  :
                            cPixBox<Dim>(aP0,aP1)
{
}

template <const int Dim> cDataGenUnTypedIm<Dim>::~cDataGenUnTypedIm()
{

}

template <const int Dim> void cDataGenUnTypedIm<Dim>::VI_VPtsSetV(const  std::vector<tPixI> & aVPt, int  aV)
{
   for (const auto & aPix : aVPt)
       VI_SetV(aPix,aV);
}

template <const int Dim> void cDataGenUnTypedIm<Dim>::VD_VPtsSetV(const  std::vector<tPixI> & aVPt,tREAL8  aV)
{
   for (const auto & aPix : aVPt)
       VI_SetV(aPix,aV);
}


template <class Type,int Dim>
cDataGenUnTypedIm<Dim> * Tpl_ReadAllocImGen(const cTplBox<int,Dim>& aBox, const cDataFileIm2D *aDFI, const cPtxd<int,Dim>& aP0)
{
    if constexpr(Dim == 2) {
        auto aDIm  =  new  cDataIm2D<Type>(aBox.P0(),aBox.P1());
        if (aDFI)
            aDIm->Read(*aDFI,aP0);
        return aDIm;
    } else if constexpr (Dim ==1) {
        auto aDIm  =  new  cDataIm1D<Type>(aBox.P0(),aBox.P1());
        return aDIm;
    } else if constexpr (Dim ==3) {
        auto aDIm  =  new  cDataIm3D<Type>(aBox.P0(),aBox.P1());
        return aDIm;
    } else {
        auto aDIm  =  new  cDataTypedIm<Type,Dim>(aBox.P0(),aBox.P1());
        return aDIm;
    }
}

template<int Dim>
static cDataGenUnTypedIm<Dim> * ReadAllocImGen(eTyNums aType, const cTplBox<int,Dim>& aBox, const cDataFileIm2D *aDFI=nullptr, const cPtxd<int,Dim>& aP0 = cPtxd<int,Dim>{})
{
    switch (aType)
    {
    case eTyNums::eTN_U_INT1 : return Tpl_ReadAllocImGen<tU_INT1>(aBox,aDFI,aP0);
    case eTyNums::eTN_U_INT2 : return Tpl_ReadAllocImGen<tU_INT2>(aBox,aDFI,aP0);
    case eTyNums::eTN_INT1   : return Tpl_ReadAllocImGen<tINT1>(aBox,aDFI,aP0);
    case eTyNums::eTN_INT2   : return Tpl_ReadAllocImGen<tINT2>(aBox,aDFI,aP0);
    case eTyNums::eTN_INT4   : return Tpl_ReadAllocImGen<tINT4>(aBox,aDFI,aP0);
    case eTyNums::eTN_REAL4  : return Tpl_ReadAllocImGen<tREAL4>(aBox,aDFI,aP0);
    default : break;
    }
    MMVII_INTERNAL_ERROR("Unhandled type in ReadIm2DGen");
    return nullptr;
}

template<int Dim>
cDataGenUnTypedIm<Dim> * AllocImGen(const cPtxd<int,Dim>& aSz, eTyNums aType)
{
    return ReadAllocImGen(aType, cTplBox<int,Dim>(aSz));
}

template<int Dim>
cDataGenUnTypedIm<Dim> * AllocImGen(const cPtxd<int,Dim>& aP0, const cPtxd<int,Dim>& aP1, eTyNums aType)
{
    return ReadAllocImGen(aType, cTplBox<int,Dim>(aP0,aP1));
}

// Instanciate previous functions
template cDataGenUnTypedIm<1> * AllocImGen(const cPtxd<int,1>& aSz, eTyNums aType);
template cDataGenUnTypedIm<2> * AllocImGen(const cPtxd<int,2>& aSz, eTyNums aType);
template cDataGenUnTypedIm<3> * AllocImGen(const cPtxd<int,3>& aSz, eTyNums aType);
template cDataGenUnTypedIm<1> * AllocImGen(const cPtxd<int,1>& aP0, const cPtxd<int,1>& aP1, eTyNums aType);
template cDataGenUnTypedIm<2> * AllocImGen(const cPtxd<int,2>& aP0, const cPtxd<int,2>& aP1, eTyNums aType);
template cDataGenUnTypedIm<3> * AllocImGen(const cPtxd<int,3>& aP0, const cPtxd<int,3>& aP1, eTyNums aType);

static cDataGenUnTypedIm<2> * ReadIm2DGen(const std::string& aName, std::optional<eTyNums> aOptType, const cBox2di  &aBox)
{
    cDataFileIm2D  aDFI = cDataFileIm2D::Create(aName,eForceGray::Yes);
    auto aSz = aDFI.Sz();
    auto aP0 = cPt2di{};
    if (! aBox.IsEmpty())
    {
        aSz = aBox.Sz();
        aP0 = aBox.P0();
    }
    auto aType = aOptType ? *aOptType : aDFI.Type();
    return ReadAllocImGen(aType,cBox2di(aSz),&aDFI,aP0);
}


cDataGenUnTypedIm<2> * ReadIm2DGen(const std::string &aName,const cBox2di& aBox)
{
    return ReadIm2DGen(aName, std::nullopt, aBox);
}


cDataGenUnTypedIm<2> * ReadIm2DGen(const std::string &aName, eTyNums aType, const cBox2di& aBox)
{
    return ReadIm2DGen(aName,std::optional(aType),aBox);
}



template <const int Dim>
double cDataGenUnTypedIm<Dim>::ClipedGetValueInterpol(const cInterpolator1D &,const cPtxd<double,Dim> &,double ,bool* ) const
{
    MMVII_INTERNAL_ERROR("Unhandled type in cDataGenUnTypedIm::GetValueAndGradInterpol");
    return 0;
}

template <const int Dim>
std::pair<tREAL8,cPtxd<double,Dim>> cDataGenUnTypedIm<Dim>::GetValueAndGradInterpol(const cDiffInterpolator1D &, const cPtxd<double, Dim> &) const
{
    MMVII_INTERNAL_ERROR("Unhandled type in cDataGenUnTypedIm::GetValueAndGradInterpol");
    return {0.,cPtxd<double,Dim>{}};
}


template <const int Dim>
cDataGenUnTypedIm<Dim>* cDataGenUnTypedIm<Dim>::AllocReSampleGen(
        const cInterpolator1D &anInterpol,
        const cDataInvertibleMapping<tREAL8, Dim> &aMap,
        const cPixBox<Dim> aBox,
        double aDefValOut) const
{
    auto aResult = AllocImGen(aBox.Sz(),this->TypeVal());

    for (auto & aPixOut : *aResult)
    {
        auto aPixIn = aMap.Inverse(MMVII::ToR(aPixOut+aBox.P0()));
        auto val = this->Inside(MMVII::ToI(aPixIn)) ? this->ClipedGetValueInterpol(anInterpol,aPixIn) : aDefValOut;
        aResult->VD_SetV(aPixOut,val);
    }
    return aResult;
}

template <const int Dim>
std::pair<cPtxd<int,Dim>,cDataGenUnTypedIm<Dim>*> cDataGenUnTypedIm<Dim>::AllocReSampleGen(
    const cInterpolator1D &anInterpol,
    const cDataInvertibleMapping<tREAL8, Dim> &aMap,
    double aDefValOut) const
{
    auto aBoxOut = aMap.BoxOfFrontier(this->ToR(),1.0).ToI();
    return {aBoxOut.P0(),this->AllocReSampleGen(anInterpol,aMap,aBoxOut,aDefValOut)};
}



/* ========================== */
/*          ::                */
/* ========================== */
    //  To test Error_Handler mecanism

static std::string MesNegSz="Negative size in rect object";
static std::string  TestErHandler;
static void TestBenchRectObjError(const std::string & aType,const std::string &  aMes,const char * aFile,int aLine)
{
   TestErHandler = aMes;
}



/* ========================== */
/*          cDataTypedIm      */
/* ========================== */


template <class Type,const int Dim>
    cDataTypedIm<Type,Dim>::cDataTypedIm(const cPtxd<int,Dim> & aP0,const cPtxd<int,Dim> & aP1,Type *aRawDataLin,eModeInitImage aModeInit) :
        cDataGenUnTypedIm<Dim>(aP0,aP1),
        mDoAlloc (aRawDataLin==0),
        mRawDataLin (mDoAlloc ? cMemManager::Alloc<Type>(NbElem())  : aRawDataLin),
        mNbElemMax  (NbElem())
{
   Init(aModeInit);
}

template <>   cDataTypedIm<tREAL8,1> * cDataTypedIm<tREAL8,1>::AllocIm(const tPix& aPix)
{
   return new cDataIm1D<tREAL8>(cPt1di(0),aPix);
}
template <>   cDataTypedIm<tREAL8,2> * cDataTypedIm<tREAL8,2>::AllocIm(const tPix& aPix)
{
   return new cDataIm2D<tREAL8>(cPt2di(0,0),aPix);
}
template <>   cDataTypedIm<tREAL8,3> * cDataTypedIm<tREAL8,3>::AllocIm(const tPix& aPix)
{
   return new cDataIm3D<tREAL8>(aPix);
}



template <class Type,const int Dim>
    void cDataTypedIm<Type,Dim>::Resize(const cPtxd<int,Dim> & aP0,const cPtxd<int,Dim> & aP1,eModeInitImage aModeInit)
{
    //  WARNING : this work because cDataGenUnTypedIm only calls cRectObj
    //     DO NOT WORK all stuff like :  this->cDataGenUnTypedIm<Dim>::cDataGenUnTypedIm(aP0,aP1);
    // static_cast<cRectObj<Dim>&>(*this) = cRectObj<Dim>(aP0,aP1);

    // this-> cPixBox<Dim>::cRectObj(aP0,aP1);

    // To call the copy constructor of cPixBox, we use a placemennt new
    // Not the best C++, but I don't success to do it other way as constructor cannot  be called explicitely
    new (static_cast<cPixBox<Dim>*>(this)) cPixBox<Dim>(aP0,aP1);

    if (cMemManager::Resize(mRawDataLin,0,mNbElemMax,0,NbElem()))
    {
        mDoAlloc = true;
    }
    Init(aModeInit);
}

template<class Type, const int Dim>
void cDataTypedIm<Type, Dim>::ToFile(const std::string &, const std::vector<std::string> &) const
{
    MMVII_INTERNAL_ERROR("Unimpemented method ToFile() (Type: " + ToStr(this->TypeVal()) + " Dim:" + std::to_string(Dim) );
}



template <class Type,const int Dim>
    cDataTypedIm<Type,Dim>::~cDataTypedIm()
{
   if (mDoAlloc)
      cMemManager::Free(mRawDataLin);
}


template <class Type,const int Dim>
        double cDataTypedIm<Type,Dim>::L1Dist(const cDataTypedIm<Type,Dim> & aI2,bool isAvg) const
{
    tPB::AssertSameArea(aI2);
    double aRes = 0.0;
    for (int aK=0 ; aK<NbElem() ; aK++)
       aRes += std::fabs(mRawDataLin[aK]-aI2.mRawDataLin[aK]);

   return isAvg ? aRes/NbElem() : aRes;
}

template <class Type,const int Dim>
        double cDataTypedIm<Type,Dim>::SqL2Dist(const cDataTypedIm<Type,Dim> & aI2,bool isAvg) const
{
    tPB::AssertSameArea(aI2);
    double aRes = 0.0;
    for (int aK=0 ; aK<NbElem() ; aK++)
       aRes += R8Square(mRawDataLin[aK]-aI2.mRawDataLin[aK]);

   return isAvg ? aRes/NbElem() : aRes;
}

template <class Type,const int Dim>
        double cDataTypedIm<Type,Dim>::L2Dist(const cDataTypedIm<Type,Dim> & aI2,bool isAvg) const
{
   return sqrt(SqL2Dist(aI2,isAvg));
}


template <class Type,const int Dim>
        double cDataTypedIm<Type,Dim>::SafeMaxRelDif(const cDataTypedIm<Type,Dim> & aI2,tREAL8 aEps) const
{
    tPB::AssertSameArea(aI2);
    double aRes = 0;
    for (int aK=0 ; aK<NbElem() ; aK++)
    {
       tREAL8 aV1 = mRawDataLin[aK];
       tREAL8 aV2 = aI2.mRawDataLin[aK];
       UpdateMax(aRes,fabs(aV1-aV2)/std::max(aEps,std::max(std::fabs(aV1),std::fabs(aV2))));
    }

   return aRes;
}

template <class Type,const int Dim>
        double cDataTypedIm<Type,Dim>::LInfDist(const cDataTypedIm<Type,Dim> & aI2) const
{
    tPB::AssertSameArea(aI2);
    double aRes = 0.0;
    for (int aK=0 ; aK<NbElem() ; aK++)
       aRes = std::max(aRes,(double)std::fabs(mRawDataLin[aK]-aI2.mRawDataLin[aK]));

   return aRes;
}



template <class Type,const int Dim>
        double cDataTypedIm<Type,Dim>::L1Norm(bool isAvg) const
{
    double aRes = 0.0;
    for (int aK=0 ; aK<NbElem() ; aK++)
       aRes += std::fabs(mRawDataLin[aK]);

   return isAvg ? aRes/NbElem() : aRes;
}
template <class Type,const int Dim>
        double cDataTypedIm<Type,Dim>::SqL2Norm(bool isAvg) const
{
    double aRes = 0.0;
    for (int aK=0 ; aK<NbElem() ; aK++)
       aRes += R8Square(mRawDataLin[aK]);

   return isAvg ? aRes/NbElem() : aRes;
}

template <class Type,const int Dim>
        double cDataTypedIm<Type,Dim>::L2Norm(bool isAvg) const
{
   return sqrt(SqL2Norm(isAvg));
}

template <class Type,const int Dim>
        double cDataTypedIm<Type,Dim>::LInfNorm() const
{
    double aRes = 0.0;
    for (int aK=0 ; aK<NbElem() ; aK++)
       aRes = std::max(aRes,(double) std::fabs(mRawDataLin[aK]));

   return aRes;
}

template <class Type,const int Dim>  Type     cDataTypedIm<Type,Dim>::MinVal() const
{
    Type aRes = mRawDataLin[0];
    for (int aK=1 ; aK<NbElem() ; aK++)
        UpdateMin(aRes,mRawDataLin[aK]);
   return aRes;
}
template <class Type,const int Dim>  Type     cDataTypedIm<Type,Dim>::MaxVal() const
{
    Type aRes = mRawDataLin[0];
    for (int aK=1 ; aK<NbElem() ; aK++)
        UpdateMax(aRes,mRawDataLin[aK]);
   return aRes;
}
template <class Type,const int Dim>  tREAL16     cDataTypedIm<Type,Dim>::SomVal() const
{
    tREAL16 aRes = mRawDataLin[0];
    for (int aK=1 ; aK<NbElem() ; aK++)
        aRes += mRawDataLin[aK];
   return aRes;
}
template <class Type,const int Dim>  tREAL16     cDataTypedIm<Type,Dim>::MoyVal() const
{
   return SomVal() / NbElem();
}



template <class Type,const int Dim> void  cDataTypedIm<Type,Dim>::DupIn(cDataTypedIm<Type,Dim> & aIm) const
{
    tPB::AssertSameSz(aIm);
    MemCopy(aIm.RawDataLin(),RawDataLin(),NbElem());
    // MMVII_INTERNAL_ASSERT_strong(mSz[aK]>=0,"");
}

template <class Type,const int Dim> void  cDataTypedIm<Type,Dim>::DupInVect(std::vector<Type> & aVec) const
{
    aVec.resize(NbElem());
    MemCopy(aVec.data(),RawDataLin(),NbElem());
    // MMVII_INTERNAL_ASSERT_strong(mSz[aK]>=0,"");
}


template <class Type,const int Dim> void  cDataTypedIm<Type,Dim>::InitCste(const Type & aVal)
{
   if (aVal==0)
   {
      InitNull();
   }
   else
   {
      // StdOut() << "xxInitCsteInitCste " << NbElem() << " " << mRawDataLin << "\n";
      for (tINT8 aK=0 ; aK< NbElem() ; aK++)
           mRawDataLin[aK] = aVal;
   }
}

template <class Type,const int Dim> void  cDataTypedIm<Type,Dim>::InitBorder(const Type & aVal)
{
   int aLarg = 1;
   if (MinAbsCoord(tPB::Sz()) > (2*aLarg))
   {
      cBorderPixBox<Dim> aBorder(this->RO(),aLarg);

      for (const auto & aP : aBorder)
      {
          mRawDataLin[tPB::IndexeLinear(aP)] = aVal;
      }
   }
   else
   {
      InitCste(aVal);
   }
}

template <class Type,const int Dim> void  cDataTypedIm<Type,Dim>::InitInteriorAndBorder(const Type & aVInt,const Type & aVB)
{
   InitCste(aVInt);
   InitBorder(aVB);
}


template <class Type,const int Dim> void  cDataTypedIm<Type,Dim>::InitRandom()
{
   for (tINT8 aK=0 ; aK< NbElem() ; aK++)
       mRawDataLin[aK] = tTraits::RandomValue();
}

template <class Type,const int Dim> void  cDataTypedIm<Type,Dim>::InitRandom(const Type & aV0,const Type &aV1)
{
   for (tINT8 aK=0 ; aK< NbElem() ; aK++)
   {
       mRawDataLin[aK] = Type(aV0 + (aV1-aV0) *RandUnif_0_1());
       if (mRawDataLin[aK]==aV1)
           mRawDataLin[aK]--;
   }
}




template <class Type,const int Dim> void  cDataTypedIm<Type,Dim>::InitRandomCenter()
{
   for (tINT8 aK=0 ; aK< NbElem() ; aK++)
       mRawDataLin[aK] = tTraits::RandomValueCenter();
}

template <class Type,const int Dim> void  cDataTypedIm<Type,Dim>::InitDirac(const cPtxd<int,Dim> & aP,const Type &  aVal)
{
    InitNull();
    mRawDataLin[tPB::IndexeLinear(aP)] = aVal;
}

template <class Type,const int Dim> void  cDataTypedIm<Type,Dim>::InitDirac(const Type &  aVal)
{
    InitDirac((tPB::mP0+tPB::mP1)/2,aVal);
}



template <class Type,const int Dim> void  cDataTypedIm<Type,Dim>::InitNull()
{
    MEM_RAZ(mRawDataLin,NbElem());
}

template <class Type,const int Dim> void  cDataTypedIm<Type,Dim>::InitId()
{
   // Check it is a square matrix
   MMVII_INTERNAL_ASSERT_bench((Dim==2)  ,"Init Id : dim !=2");
   for (int aK=0 ; aK<Dim ; aK++)
   {
       MMVII_INTERNAL_ASSERT_bench((P0()[aK]==0)  ,"Init Id P0!= (0,0)");
   }
   MMVII_INTERNAL_ASSERT_bench((P1()[0]==P1()[1])  ,"Init Id, non square image");

   InitNull();
   for (int aK=0 ;  aK<NbElem()  ; aK += P1()[0]+1)
   {
       mRawDataLin[aK] = 1;
   }
}

template <class Type,const int Dim> void  cDataTypedIm<Type,Dim>::Init(eModeInitImage aMode)
{
    switch(aMode)
    {
       case eModeInitImage::eMIA_Rand        : InitRandom(); return;
       case eModeInitImage::eMIA_RandCenter  : InitRandomCenter(); return;
       case eModeInitImage::eMIA_Null        : InitNull(); return;
       case eModeInitImage::eMIA_V1          : InitCste(1); return;
       case eModeInitImage::eMIA_MatrixId    : InitId(); return;
       case eModeInitImage::eMIA_NoInit      : ;
    }
}


template <class Type,const int Dim> int  cDataTypedIm<Type,Dim>::VI_GetV(const cPtxd<int,Dim> & aP)  const
{
    tPB::AssertInside(aP);
    return round_ni(mRawDataLin[tPB::IndexeLinear(aP)]);
}
template <class Type,const int Dim> double  cDataTypedIm<Type,Dim>::VD_GetV(const cPtxd<int,Dim> & aP)  const
{
    tPB::AssertInside(aP);
    return mRawDataLin[tPB::IndexeLinear(aP)];
}
template <class Type,const int Dim> void  cDataTypedIm<Type,Dim>::VI_SetV(const cPtxd<int,Dim> & aP,const int & aV)
{
    tPB::AssertInside(aP);
    mRawDataLin[tPB::IndexeLinear(aP)] = tNumTrait<Type>::Trunc(aV);
}

template <class Type,const int Dim> void  cDataTypedIm<Type,Dim>::VD_SetV(const cPtxd<int,Dim> & aP,const double & aV)
{
    tPB::AssertInside(aP);
    MMVII_INTERNAL_ASSERT_tiny(tNumTrait<Type>::ValueOk(aV),"Bad Value in VD_SetV");
    mRawDataLin[tPB::IndexeLinear(aP)] = tNumTrait<Type>::RoundNearestToType(aV);
}

template<class Type,const int Dim> void  cDataTypedIm<Type,Dim>::ChSignIn(cDataTypedIm<Type,Dim> & aRes) const
{
   this->AssertSameArea(aRes);
   auto  aOut =  aRes.mRawDataLin;
   auto  aIn =   mRawDataLin;
   auto aNbElem = NbElem();

// msvc++ :  disable: warning C4146: unary minus operator applied to unsigned type, result still unsigned
#ifdef _WIN32
# pragma warning( push )
# pragma warning( disable : 4146 )
#endif

   for (int aX=0 ; aX<aNbElem ; aX++)
       aOut[aX] = -aIn[aX];

#ifdef _WIN32
# pragma warning( pop )
#endif
}


/*
template class cDataTypedIm<tREAL4,1>;
template class cDataTypedIm<tREAL4,2>;
template class cDataTypedIm<tREAL4,3>;
*/

template class cDataGenUnTypedIm<1>;
template class cDataGenUnTypedIm<2>;
template class cDataGenUnTypedIm<3>;
template class cDataGenUnTypedIm<4>;
template class cDataGenUnTypedIm<5>;





#define MACRO_INSTANTIATE_cDataTypedIm(aType)\
template class cDataTypedIm<aType,1>;\
template class cDataTypedIm<aType,2>;\
template class cDataTypedIm<aType,3>;\
template class cDataTypedIm<aType,4>;\
template class cDataTypedIm<aType,5>;


MACRO_INSTANTIATE_cDataTypedIm(tINT1)
MACRO_INSTANTIATE_cDataTypedIm(tINT2)
MACRO_INSTANTIATE_cDataTypedIm(tINT4)
MACRO_INSTANTIATE_cDataTypedIm(tINT8)

MACRO_INSTANTIATE_cDataTypedIm(tU_INT1)
MACRO_INSTANTIATE_cDataTypedIm(tU_INT2)
MACRO_INSTANTIATE_cDataTypedIm(tU_INT4)

MACRO_INSTANTIATE_cDataTypedIm(tREAL4)
MACRO_INSTANTIATE_cDataTypedIm(tREAL8)
MACRO_INSTANTIATE_cDataTypedIm(tREAL16)
/*
*/


/* ========================== */
/*          cBenchBaseImage   */
/* ========================== */


    ///  cBenchBaseImage
class     cBenchBaseImage
{
    public :
    static void DoBenchBI();
    static void DoBenchRO();
};

void cBenchBaseImage::DoBenchBI()
{
    cMemState  aState = cMemManager::CurState() ;
    {

        // cDataTypedIm<tREAL4,2> aBI(cPt2di(2,3),cPt2di(10,9));
    }
    cMemManager::CheckRestoration(aState);
}


void cBenchBaseImage::DoBenchRO()
{
   {
      cPixBox<1>  aR(cPt1di(2),cPt1di(10));
      MMVII_INTERNAL_ASSERT_bench(cPt1di(aR.Sz()) ==cPt1di(8),"Bench sz RectObj");
      MMVII_INTERNAL_ASSERT_bench(aR.NbElem()==8,"Bench sz RectObj");
   }

   // Test the SetErrorHandler mecanism, abality to recover on error
   //     => do it only if verification is high;  else the error will not be detected and check cannot work
   if (The_MMVII_DebugLevel >= The_MMVII_DebugLevel_InternalError_tiny)
   {
       MMVII_SetErrorHandler(TestBenchRectObjError);
       cPixBox<1>  aR(cPt1di(10),cPt1di(0));
       MMVII_RestoreDefaultHandle();
       MMVII_INTERNAL_ASSERT_bench(MesNegSz==TestErHandler,"Handler mechanism");
   }
   {
      cPixBox<2>  aR(cPt2di(2,3),cPt2di(10,9));
      MMVII_INTERNAL_ASSERT_bench(cPt2di(aR.Sz()) ==cPt2di(8,6),"Bench sz RectObj");
      MMVII_INTERNAL_ASSERT_bench(aR.NbElem()==48,"Bench sz RectObj");

      MMVII_INTERNAL_ASSERT_bench(aR.Inside(cPt2di(2,8)),"Bench inside rect");
      MMVII_INTERNAL_ASSERT_bench(aR.Inside(cPt2di(9,3)),"Bench inside rect");
      MMVII_INTERNAL_ASSERT_bench(!aR.Inside(cPt2di(2,9)),"Bench inside rect");
      MMVII_INTERNAL_ASSERT_bench(!aR.Inside(cPt2di(1,8)),"Bench inside rect");
      MMVII_INTERNAL_ASSERT_bench(!aR.Inside(cPt2di(10,3)),"Bench inside rect");
      MMVII_INTERNAL_ASSERT_bench(!aR.Inside(cPt2di(9,2)),"Bench inside rect");
      MMVII_INTERNAL_ASSERT_bench(!aR.Inside(cPt2di(10,2)),"Bench inside rect");
   }
}

    //----------  ::  external call ---------

void BenchBaseImage()
{
    cBenchBaseImage::DoBenchBI();
}

void BenchRectObj()
{
    cBenchBaseImage::DoBenchRO();
}




};
