#ifndef  _MMVII_Tpl_NormImages_H_
#define  _MMVII_Tpl_NormImages_H_

#include "MMVII_Matrix.h"

namespace MMVII
{



         /*****************************************************/
         /*                                                   */
         /*          NORMALIZATION                            */
         /*                                                   */
         /*****************************************************/

         ///  Compute Average and Standard Deviation of 1 image
template<class T,int Dim>
cComputeStdDev<tREAL8>  ComputedAvgDev(const cDataTypedIm<T,Dim> & aIm,tREAL8 aEpsilon)
{
     cComputeStdDev<tREAL8>  aCSD;
     for (int aK=0 ; aK<aIm.NbElem() ; aK++)
     {
         aCSD.Add(1.0,aIm.GetRDL(aK));
     }
     aCSD = aCSD.Normalize(aEpsilon);

     return aCSD;
}

/// In place, replace value by normalized value using Avg & StdDev
template<class T,int Dim>
void GenNormalizedAvgDev(cDataTypedIm<T,Dim> & aIm,tREAL8 aEpsilon,tREAL8 aMul=1.0)
{


    cComputeStdDev<tREAL8>  aCSD = ComputedAvgDev(aIm,aEpsilon);

    for (int aK=0 ; aK<aIm.NbElem() ; aK++)
    {
        T& aV = aIm.GetRDL(aK);
        aV = aCSD.NormalizedVal(aV) * aMul;
    }
}



template<class T>  cIm2D<T> Im2DNormalizedAvgDev(const cIm2D<T> & aIm,tREAL8 aEpsilon=1e-5)
{
    cIm2D<T> aRes = aIm.Dup();
    GenNormalizedAvgDev(aRes.DIm(),aEpsilon);
    return aRes;
}

         /*****************************************************/
         /*                                                   */
         /*         8-Bits transo, usefull for visu           */
         /*                                                   */
         /*****************************************************/


    /* ======================================================================================= */
    /* ==============================  GAUSSIAN NORMALIZATION ================================ */
    /* ======================================================================================= */


 /** Tranformat Im in 8Bit Res, using a gaussian modelization : compute avg & ast-dev,
  *  then use error function to map in :
  *
  *          [-1,1] => [TargetMin,TargetMax]
  */

template<class TOut,class TIn,int Dim>
       void FixDynByGaussModel
            (
               cDataTypedIm<TOut,Dim> & aRes,
               tREAL8 aTargetMin,
               tREAL8 aTargetMax,
               const cDataTypedIm<TIn,Dim> & aIm,
               tREAL8 aEpsilon=1e-5
             )
{
    aRes.AssertSameArea(aIm);

    cComputeStdDev<tREAL8>  aCSD = ComputedAvgDev(aIm,aEpsilon);

    tREAL8 aVCenter = (aTargetMin+aTargetMax) / 2.0;
    tREAL8 aDyn = (aTargetMax-aTargetMin) / 2.0;

    for (int aK=0 ; aK<aIm.NbElem() ; aK++)
    {
        tREAL8 aV = aIm.GetRDL(aK);
        aV = aCSD.NormalizedVal(aV);
        aV = aVCenter + aDyn * std::erf(aV);
        aRes.GetRDL(aK) = tNumTrait<TOut>::Trunc(aV);
    }
}

template<class T,int Dim>
      void Make8BByGaussModel(cDataTypedIm<tU_INT1,Dim> & aRes,const cDataTypedIm<T,Dim> & aIm,tREAL8 aEpsilon=1e-5)
{
     FixDynByGaussModel(aRes,0,256,aIm,aEpsilon);
}


template<class T> cIm2D<tU_INT1> Im2Make8BByGaussModel(const cDataIm2D<T> & aIm,tREAL8 aEpsilon=1e-5)
{
     cIm2D<tU_INT1>  aRes(aIm.Sz());
     Make8BByGaussModel(aRes.DIm(),aIm,aEpsilon);
     return aRes;
}


      /* ======================================================================================= */
      /* ============================== MIN-MAX NORMALIZATION ================================== */
      /* ======================================================================================= */


/** Tranformat Im in 8Bit Res, using a linear transfo :
 *
 *          [Min,Max] => [TargetMin,TargetMax]
*/

template<class TOut,class TIn,int Dim>
     void FixDynBByMinMax
         (
             cDataTypedIm<TOut,Dim> & aRes,
             const cDataTypedIm<TIn,Dim> & aIm,
             tREAL8  aTargetMin,
             tREAL8  aTargetMax,
             tREAL8  aEpsilon=1e-5
        )
{
  aRes.AssertSameArea(aIm);

  TIn aVMin,aVMax;
  GetBounds(aVMin,aVMax,aIm);

  tREAL8 aMul = (aTargetMax-aTargetMin) / std::max(aEpsilon,tREAL8(aVMax-aVMin));
  for (int aK=0 ; aK<aIm.NbElem() ; aK++)
  {
      tREAL8 aV = aIm.GetRDL(aK);
      aV =  aTargetMin+ (aV-aVMin) *  aMul  ;
      aRes.GetRDL(aK) = tNumTrait<TOut>::Trunc(aV);
  }
}

template<class T,int Dim>
          void Make8BByMinMax
              (
                  cDataTypedIm<tU_INT1,Dim> & aRes,
                  const cDataTypedIm<T,Dim> & aIm,
                  tREAL8 aEpsilon=1e-5
             )
{
   FixDynBByMinMax(aRes,aIm,0,256,aEpsilon);
}

template<class T> cIm2D<tU_INT1> Im2Make8BByMinMax(const cDataIm2D<T> & aIm,tREAL8 aEpsilon=1e-5)
{
     cIm2D<tU_INT1>  aRes(aIm.Sz());
     Make8BByMinMax(aRes.DIm(),aIm,aEpsilon);
     return aRes;
}

          /* ======================================================================================= */
          /* ========== MAX NORMALIZATION (only mult, may be pertinent for >=0 value) ============== */
          /* ======================================================================================= */

          /** Tranformat Im in 8Bit Res, using a linear transfo :
           *
           *          [Min,Max] => [TargetMin,TargetMax]
          */

template<class TOut,class TIn,int Dim>
void FixDynBByMax
     (
          cDataTypedIm<TOut,Dim> & aRes,
          const cDataTypedIm<TIn,Dim> & aIm,
          tREAL8  aTargetMax,
          tREAL8  aEpsilon=1e-5
    )
{
      aRes.AssertSameArea(aIm);

      TIn aUnusedMin,aVMax;
      GetBounds(aUnusedMin,aVMax,aIm);

      tREAL8 aMul = (aTargetMax) / std::max(aEpsilon,tREAL8(aVMax));
      for (int aK=0 ; aK<aIm.NbElem() ; aK++)
      {
          tREAL8 aV = aIm.GetRDL(aK);
          aV = aV *  aMul  ;
          aRes.GetRDL(aK) = tNumTrait<TOut>::Trunc(aV);
      }
}


template<class T>
          cIm1D<T> Im1DFixDynByMax
              (
                  const cDataIm1D<T> & aIm,
                  tREAL8  aTargetMax,
                  tREAL8 aEpsilon=1e-5
             )
{
    cIm1D<T>  aRes(aIm.Sz());
    FixDynBByMax(aRes.DIm(),aIm,aTargetMax,aEpsilon);
    return aRes;
}


     /* ======================================================================================= */
     /* ========================= HISTOGRAM-EQUALISATIO? NORMALIZATION ======================== */
     /* ======================================================================================= */

 /** Tranformat Im in 8Bit Res, using a gaussian modelization : compute avg & ast-dev,
 *  then use error function to map in [-1,1] -> [0,255]
 */

 template<class T,int Dim>
      void Make8BByHistoEq(cDataTypedIm<tU_INT1,Dim> & aRes,const cDataTypedIm<T,Dim> & aIm,int aNbVal=100000,tREAL8 aEpsilon=1e-5)
 {
   aRes.AssertSameArea(aIm);

   T aVMin,aVMax;
   GetBounds(aVMin,aVMax,aIm);
   tREAL8 aMul = aNbVal / std::max(aEpsilon,tREAL8(aVMax-aVMin));

   cHistoCumul<tREAL8,tREAL8> aHC(aNbVal);

   for (int aK=0 ; aK<aIm.NbElem() ; aK++)
   {
       tREAL8 aV = aIm.GetRDL(aK);
       aV = (aV-aVMin) *  aMul  ;
       aHC.AddV(round_ni(aV),1.0);
   }
   aHC.MakeCumul();

   for (int aK=0 ; aK<aIm.NbElem() ; aK++)
   {
       tREAL8 aV = aIm.GetRDL(aK);
       aV = (aV-aVMin) *  aMul  ;
       aV = 255.0 * aHC.PropCumul(aV);
       aRes.GetRDL(aK) = tNumTrait<tU_INT1>::Trunc(aV);
   }
 }

template<class T> cIm2D<tU_INT1> Im2Make8BByHistoEq(const cDataIm2D<T> & aIm,int aNbLine=100000,tREAL8 aEpsilon=1e-5)
{
   cIm2D<tU_INT1>  aRes(aIm.Sz());
   Make8BByHistoEq(aRes.DIm(),aIm,aNbLine,aEpsilon);
   return aRes;
}



};

#endif  //  _MMVII_Tpl_NormImages_H_
