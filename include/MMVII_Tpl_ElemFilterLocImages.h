#ifndef  _MMVII_Tpl_ElemFilterLocImages_
#define  _MMVII_Tpl_ElemFilterLocImages_

#include "MMVII_Matrix.h"



/** \file   MMVII_Tpl_ElemFilterLocImages.h
    \brief  Implemantation of "convulation" operators on images,
            they are "elementary" because there is no "clever" optimization
            like can be done with fast average, fast max ....
*/



namespace MMVII
{

/// Interface for vector to aggregate vals
/*
template <class Type>    class cVectAddVal
{
    public :
       std::vector<Type>  mVVals;
       void AddVal(Type aVal) {mVVals.push_back(aVal);}
};

template<class TypeVal,class TypeAgreg >
    void  PushValInAgreg (TypeAgreg & anAgreg,const cDataIm2D<Type> & aI2,const cBox2di & aBoxFul)
{
    cBox2di aBoc
    for (const auto &)
}
        */

template<class Type>   Type  KthValLoc (const cDataIm2D<Type> & aIm,const cPt2di & aP0, const cRect2 & aNeigh,tREAL8 aProp,Type aDef)
{
    static   std::vector<tREAL8> aVVals;
    aVVals.clear();

    for (const auto & aPN : aNeigh)
    {
        cPt2di aPt = aP0 + aPN;
        if (aIm.Inside(aPt))
        {
            aVVals.push_back(aIm.GetV(aPt));
        }
    }
    if (aVVals.empty())
        return aDef;
    return NC_KthVal(aVVals,aProp);
}

};

#endif  //  _MMVII_Tpl_ElemFilterLocImages_
