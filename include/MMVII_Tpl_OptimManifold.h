#include "MMVII_HeuristikOpt.h"
#include "MMVII_Mappings.h"

namespace MMVII
{

/** @file
 *
 *    Contains class for optimizing, on a bounded manifold, "any" function regular enough.
 *
 *    The typicall manifold can be : sphere, rotation, torus, space of poses ....
 *
 *    The function has to +- regular for the approach to work, but not necessarily differentiable nor
 *  modelized as a least square.
 *
 *    The approach is :
 *        - sample the manifold
 *        - extract a certain number of local mimimal on these samples
 *        - refine these minimal on the "tangent space" using the cOptimByStep-strategy,
 *        because the neigboorhoud can be assimilated to a vectorial space
 *
 *   =======================================================================================================
 *
 *    The optimization is done in the template class cTplOptDisc_OnManifold<tParamManif>. The class
 *  tParamManif describe the manifold and how to exlporate it. There is  several requirement.
 *
 *    It must define dimension and types :
 *
 *      - define the dimension of the tangent space "DimTgtSp=DT"
 *      - define the point tangent space by  :      typedef cPtxd<tREAL8,DimTgtSp> tPtTgt;
 *        (not mandatory but recommanded to uniformize the manipulation)
 *      - define the type of object stored in the manifold by the typedef "tObj"
 *      - define  the type allowing to compute the tangent space by  "tDescrTgtSp"
 *
 *   It must provide the following method to sample the manifold :
 *
 *      - int NbObjDisc() const : return the number of samples used on the manifold
 *      - tObj KthObjDisc(size_t aK) const  : return the kTh object
 *      - tREAL8 Density() const : return an estimation of the density of the sampling (not accurate)
 *      - tREAL8 SqDist(const tObj & aP1,const tObj & aP2) const : return the square distance between 2
 *         object of the manifold
 *
 *  It must provide the following method to describe the tangent space :
 *
 *     - tDescrTgtSp * CreateTgtSpace(const tObj & aPt) const , return the object for manipulating tangent space in
 *       point Pt
 *     - void DeleteObjFromTgtSp(tDescrTgtSp * anObj) const, deallocate the tangent space
 *
 *     - tObj GetObjFromTgSpace(const tDescrTgtSp & , const tPtTgt& aPTgt) const : given a descriptor of
 *      tangent point and point of R^DT, return the point of manifold corresponding tp aPTgt
 *
 *     =======================================================================================================
 *
 *  The classes  cSph3_OptimDisc and  cRot_OptimDisc are two example of such classes.
 *
 *
 *      "cSph3_OptimDisc" : is used to desricibe the sphere. The object used for describing the
 *      tangent space to point U is  "cP3dNormWithUK" as it contains two vectors V,W such that
 *      (U,V,W) is ortho-normal
 *
 *      "cRot_OptimDisc" : is used to describe the set of rotation. The object used to describe
 *      the tangent space to a rotation "R" is simply the "R" itsfelf because the neighboors of "R"
 *      can be computed as R*(I+dR) where "I+dR" is a small rotation.
 *
 *       =======================================================================================================
 *
 *     If "M1" and "M2" are two manifold, the cartesian product  "M1 X M2" is also a manifold. The template class
 *   "cCartProduct_OptimDisc" allow to create automatically a descritor of "M1 X M2" from a descriptor of M1 and M2.
 *   For example, a poses being made of unit vector and a rotation, the manifold of the space of poses can be
 *   constructed as this :
 *
 *            typedef cCartProduct_OptimDisc<cSph3_OptimDisc,cRot_OptimDisc> tPoseCart_OptimDisc;
 *
 *
 */


/* ******************************************** */
/*                                              */
/*              cSph3_OptimDisc                 */
/*                                              */
/* ******************************************** */

/** Class for optimizing function on the sphere */

class cSph3_OptimDisc
{
   public :
       typedef cPt3dr tObj;
       static constexpr int DimTgtSp =2;
       typedef cP3dNormWithUK  tDescrTgtSp;   
       typedef cPtxd<tREAL8,DimTgtSp> tPtTgt;


       cSph3_OptimDisc(int aNbSS,bool isSphereProj)  :
            mSampSph(aNbSS,isSphereProj)
       {
       }

       int NbObjDisc() const {return mSampSph.NbSamples();}
       tObj KthObjDisc(size_t aK) const { return mSampSph.KthPt(aK); }

       /// Distance between point can be basic or projective (P ~ -P) according to Sphere
       tREAL8 SqDist(const tObj & aP1,const tObj & aP2) const
       {
           return mSampSph.SqDist(aP1,aP2);
       }
       tREAL8 Density() const {return 1.0/NbObjDisc();}


       tDescrTgtSp * CreateTgtSpace(const tObj & aPt) const
       {
           return  new tDescrTgtSp(aPt,"","");
       }
       void DeleteObjFromTgtSp(tDescrTgtSp * anObj) const
       {
           delete anObj;
       }
       tObj GetObjFromTgSpace(const tDescrTgtSp & aP3N, const tPtTgt& aPTgt) const
       {
           return aP3N.GetPNorm(aPTgt);
       }



    private :
       cSampleSphere3D  mSampSph;
};

/* ******************************************** */
/*                                              */
/*              cRot_OptimDisc                  */
/*                                              */
/* ******************************************** */


class cRot_OptimDisc
{
    public :
       static constexpr int DimTgtSp=3;
       typedef tRotR tObj;
       typedef tRotR  tDescrTgtSp;
       typedef cPtxd<tREAL8,DimTgtSp> tPtTgt;



       cRot_OptimDisc(int aNbSQuat) :
          mSampQuat(aNbSQuat,true)
       {
       }

       int NbObjDisc() const {return mSampQuat.NbRot();}
       tRotR KthObjDisc(size_t aK) const { return mSampQuat.KthRot(aK); }

       tREAL8 SqDist(const tObj & aR1,const tObj & aR2) const
       {
           //  Div 2.0 => +- to have it homogeneous to radian for rot close 2 ID
           return aR1.Mat().DIm().SqL2Dist(aR2.Mat().DIm()) / 2.0;
       }

       tDescrTgtSp * CreateTgtSpace(const tObj & aRot) const
       {
           return new tDescrTgtSp(aRot);
       }
       void DeleteObjFromTgtSp(tDescrTgtSp * anObj) const
       {
           delete anObj;
       }

       tObj GetObjFromTgSpace(const tDescrTgtSp & aRot, const tPtTgt& aPTgt) const
       {
           return aRot * tRotR::RotFromWPK(aPTgt);
       }

       tREAL8 Density() const {return 1.0/NbObjDisc();}


    private :
       cSampleQuat  mSampQuat;
};

/* ******************************************** */
/*                                              */
/*              cCartProduct_OptimDisc          */
/*                                              */
/* ******************************************** */


template <class T1,class T2> class cCartProduct_OptimDisc
{
    public :
       static constexpr int DimTgtSp= T1::DimTgtSp + T2::DimTgtSp;
       typedef std::pair<typename T1::tObj,typename T2::tObj> tObj; 
       typedef std::pair<typename T1::tDescrTgtSp*,typename T2::tDescrTgtSp*>  tDescrTgtSp;
       typedef cPtxd<tREAL8,DimTgtSp> tPtTgt;

       cCartProduct_OptimDisc(const T1& aT1,const T2& aT2) :
           mT1 (aT1),
           mNb1 (mT1.NbObjDisc()),
           mT2 (aT2),
           mNb2 (mT2.NbObjDisc())
       {
       }

       int NbObjDisc() const {return mNb1 * mNb2;}

       tREAL8 Density() const {return (mT1.Density()+mT2.Density()) / 2.0 ;}

       tObj KthObjDisc(size_t aK) const
       {
           size_t aK1 = aK % mNb1;
           size_t aK2 = aK / mNb1;
           return tObj(mT1.KthObjDisc(aK1),mT2.KthObjDisc(aK2));

       }

       tREAL8 SqDist(const tObj & aPair1,const tObj & aPair2) const
       {
           return     mT1.SqDist(aPair1.first,aPair2.first)
                   +  mT2.SqDist(aPair1.second,aPair2.second);
       }

       tDescrTgtSp *  CreateTgtSpace(const tObj & aPair) const
       {
           return new tDescrTgtSp(mT1.CreateTgtSpace(aPair.first),mT2.CreateTgtSpace(aPair.second));
       }

       void DeleteObjFromTgtSp(tDescrTgtSp * anObj) const
       {
           if (anObj)
           {
              mT1.DeleteObjFromTgtSp(anObj->first);
              mT2.DeleteObjFromTgtSp(anObj->second);

              delete anObj;
           }
       }


       tObj GetObjFromTgSpace(const tDescrTgtSp & aPair, const tPtTgt& aPTgt) const
       {
           typename T1::tPtTgt aPTg1(aPTgt.PtRawData());
           typename T2::tPtTgt aPTg2(aPTgt.PtRawData()+T1::DimTgtSp);

           return tObj
                   (
                       mT1.GetObjFromTgSpace(*aPair.first ,aPTg1),
                       mT2.GetObjFromTgSpace(*aPair.second,aPTg2)
                   );
       }

    private :
       T1 mT1;
       size_t mNb1;
       T2 mT2;
       size_t mNb2;


};

typedef cCartProduct_OptimDisc<cSph3_OptimDisc,cRot_OptimDisc> tPoseCart_OptimDisc;

inline tPoseCart_OptimDisc  Pose_OptimDisc(int aNbTr,bool isSPhProj,int aNbRot)
{
    cSph3_OptimDisc     aSph3Desc(aNbTr,isSPhProj);
    // Object that allow to sample rot-3 and explore tangent to a given rot
    cRot_OptimDisc      aRotDesc(aNbRot);

    // Object that allow to sample pose and explore tangent to a given pose
    return  tPoseCart_OptimDisc(aSph3Desc,aRotDesc);
}
/* ******************************************** */
/*                                              */
/*              cTplOptDisc_OnManifold          */
/*                                              */
/* ******************************************** */


template <class tParamManif> class cTplOptDisc_OnManifold : public cDataMapping<tREAL8,tParamManif::DimTgtSp,1>
{
     public :

         static constexpr int DimTgtSp = tParamManif::DimTgtSp;

         typedef cPtxd<tREAL8,DimTgtSp>          tPtTgt;
         typedef typename tParamManif::tObj        tObj;
         typedef typename tParamManif::tDescrTgtSp tDescrTgtSp;
         typedef cOptDiscScorer<tObj>        tScorer;
         typedef cWhichMin<tObj,tREAL8>  tResCompute;

        //  cTplOptDisc_OnManifold(const tParamGen & , const  tScorer & aScorer);
         cTplOptDisc_OnManifold
              (
                  const tParamManif & aParam,
                  const  tScorer & aScorer
               ) :
                    mManifDesc   (aParam),
                    mScorer      (aScorer),
                    mNbOjDisc    (mManifDesc.NbObjDisc()),
                    mDescTgtSp   (nullptr),
                    mMinObj      (tObj(),1e40),
                    mThsScoreEnd (-1.0)
         {
            mVObjDisc.reserve(mNbOjDisc);
            mVScoreDisc.reserve(mNbOjDisc);
            mVObjIsFree.reserve(mNbOjDisc);

            for (size_t aKObj=0 ; aKObj<mNbOjDisc ; aKObj++)
            {
                mVObjDisc.push_back(mManifDesc.KthObjDisc(aKObj));
                mVScoreDisc.push_back(mScorer.Score(mVObjDisc.back()));
                mVObjIsFree.push_back(true);
            }
         }
         const cWhichMin<tObj,tREAL8> & GetSol() const {return mMinObj;}
         cPt1dr Value(const tPtTgt& aPt) const override
         {
             tObj anObj = mManifDesc.GetObjFromTgSpace(*mDescTgtSp,aPt);

             return cPt1dr(mScorer.Score(anObj));
         }



         //void  ComputeSol(tREAL8 aDistNeigh,tREAL8 anEpsilon,int aNbTest);
             void ComputeSol
             (
                 tREAL8 aDistNeigh,
                 tREAL8 anEpsilon,
                 int aNbTest
             )
         {
             for (int aKTest = 0 ; aKTest<aNbTest; aKTest++)
             {
                 int aKMin = ExtractMin(aDistNeigh);
                 if (aKMin<0)
                 {
                     MMVII_INTERNAL_ASSERT_always(aKTest!=0,"No cdt cTplOptDisc_OnManifold<tParamGen>::ComputeSol");
                     aKTest = aNbTest;
                 }
                 else
                 {
                     mManifDesc.DeleteObjFromTgtSp(mDescTgtSp);
                     mDescTgtSp = mManifDesc.CreateTgtSpace(mVObjDisc.at(aKMin));

                     cOptimByStep<DimTgtSp> anOptim(*this,true,1.0);
                     auto [aScOpt,aPtTgt] = anOptim.Optim(tPtTgt::PCste(0),mManifDesc.Density(),anEpsilon,1/sqrt(2.0));
                     tObj anObj = mManifDesc.GetObjFromTgSpace(*mDescTgtSp,aPtTgt);

                     if (aKTest==0)
                     {
                        mThsScoreEnd = mScorer.ComputeThreshold(aScOpt,anObj);
                     }
                     else
                     {
                         if (aScOpt>mThsScoreEnd)
                         {
                            aKTest = aNbTest;
                         }
                     }

                     mMinObj.Add(anObj,aScOpt);
                 }
             }
             mManifDesc.DeleteObjFromTgtSp(mDescTgtSp);

         }


    private :

         /// extract the min candidat and mark as not free those at distance < aDistMin
         int ExtractMin(tREAL8 aDistMin)
         {
             cWhichMin<int,tREAL8> aMinK(-1,1e20);

             // extract the best scores cdt that is still free
             for (size_t aK=0 ; aK<mNbOjDisc ; aK++)
             {
                  if (mVObjIsFree.at(aK) )
                      aMinK.Add(aK, mVScoreDisc.at(aK));

             }

             int aRes = aMinK.IndexExtre();

             // is it exist, mark the neighbooring candidate as not free
             if (aRes >=0)
             {
                 tREAL8 aD2 = Square(aDistMin);
                 tObj anObj = mVObjDisc.at(aRes);

                 for (size_t aK=0 ; aK<mNbOjDisc ; aK++)
                 {
                     if (mManifDesc.SqDist(anObj,mVObjDisc.at(aK)) < aD2)
                        mVObjIsFree.at(aK) = false;
                 }


             }

             return aRes;
         }

         void  ComputeValDisclInit();

         tParamManif             mManifDesc;
         const tScorer &       mScorer;
         size_t                mNbOjDisc;
         std::vector<tObj>     mVObjDisc;
         std::vector<tREAL8>   mVScoreDisc;
         std::vector<bool  >   mVObjIsFree;
         tDescrTgtSp *         mDescTgtSp;
         cWhichMin<tObj,tREAL8> mMinObj;
         tREAL8                 mThsScoreEnd;

};



};

