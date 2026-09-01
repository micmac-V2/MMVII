
#include "MMVII_Tpl_Images.h"

#include "MMVII_SysSurR.h"

using namespace NS_SymbolicDerivative;
using namespace MMVII;

namespace MMVII
{

/* ************************************************************ */
/*                                                              */
/*                cREAL8_RWAdapt                                */
/*                                                              */
/* ************************************************************ */

template <class Type>
  cREAL8_RWAdapt<Type>::cREAL8_RWAdapt(const cBasicWeighter<tREAL8> * aRW) :
    mRW (aRW)
{
}


template <class Type>  typename cREAL8_RWAdapt<Type>::tStdVect cREAL8_RWAdapt<Type>::WeightOfResidual(const tStdVect & aVIn) const
{
     std::vector<tREAL8> aRV;
     Convert(aRV,aVIn);

     aRV = mRW->WeightOfResidual(aRV);
     tStdVect aVOut;

     return Convert(aVOut,aRV);
}



/* ************************************************************ */
/*                                                              */
/*                   cBasicWeighter                             */
/*                                                              */
/* ************************************************************ */

template <class Type>  cBasicWeighter<Type>::cBasicWeighter(const Type & aVal) :
    mVal (aVal)
{
}

template <class Type>  cBasicWeighter<Type>::~cBasicWeighter()
{
}

template <class Type>  std::vector<Type>  cBasicWeighter<Type>::WeightOfResidual(const tStdVect & aVResidual) const
{
        return tStdVect(aVResidual.size(),mVal);
}

/* ************************************************************ */
/*                                                              */
/*                cExplicitWeighter                             */
/*                                                              */
/* ************************************************************ */

template <class Type>
cResidualWeighterExplicit<Type>::cResidualWeighterExplicit(bool isSigmas, const tStdVect & aData) :
    mSigmas{}, mWeights{}
{
    tStdVect aDataInv {};
    std::for_each
    (
         aData.begin(), aData.end(),
         [&](Type aValue) { aDataInv.push_back( isSigmas ? 1/Square(aValue) : 1/std::sqrt(aValue) ); }
    );
    if (isSigmas)
    {
        mSigmas = aData;
        mWeights = aDataInv;
    }
    else
    {
        mSigmas = aDataInv;
        mWeights = aData;
    }
}

template <class Type>
std::vector<Type> cResidualWeighterExplicit<Type>::WeightOfResidual(const std::vector<Type> &aVResidual) const
{
    MMVII_INTERNAL_ASSERT_tiny(mWeights.size() == aVResidual.size(), "Number of weights does not correpond to number of residuals");
    return mWeights;
}

template <class Type>
const std::vector<Type> & cResidualWeighterExplicit<Type>::getSigmas() const
{
    return mSigmas;
}

template <class Type>
std::vector<Type> & cResidualWeighterExplicit<Type>::getSigmas()
{
    return mSigmas;
}

template <class Type>
const std::vector<Type> & cResidualWeighterExplicit<Type>::geWeights() const
{
    return mWeights;
}

template <class Type>
std::vector<Type> & cResidualWeighterExplicit<Type>::geWeights()
{
    return mWeights;
}

template <class Type>
void cResidualWeighterExplicit<Type>::addSigma(Type aSigma)
{
    mSigmas.push_back(aSigma);
    mWeights.push_back(1/Square(aSigma));
}

template <class Type>
void cResidualWeighterExplicit<Type>::addWeight(Type aWeight)
{
    mWeights.push_back(aWeight);
    mSigmas.push_back(1/std::sqrt(aWeight));
}

#define INSTANTIATE_RESOLSYSNL(TYPE)\
template class  cREAL8_RWAdapt<TYPE>;\
template class  cBasicWeighter<TYPE>;\
template class  cResidualWeighterExplicit<TYPE>;

INSTANTIATE_RESOLSYSNL(tREAL4)
INSTANTIATE_RESOLSYSNL(tREAL8)
INSTANTIATE_RESOLSYSNL(tREAL16)


};
