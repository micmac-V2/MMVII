//#include "MMVII_PCSens.h"

// Test commit

#include "cAppliNewFranges.h"

#include "MMVII_Tpl_Graph_SubGraph.h"
#include "MMVII_Tpl_GraphAlgo_SPCC.h"
//#include "MMVII_Tpl_Images.h"


namespace MMVII
{

namespace NS_FrangesDetect
{


/* ******************************************** */
/*                                              */
/*               tGraphReg                      */
/*                                              */
/* ******************************************** */
static constexpr tREAL8 mMaxValForced = 1e5;
static constexpr tREAL8 mMaxValProhibited = 0.0;

static constexpr tREAL8 mCostGeomForced     = 0.0;
static constexpr tREAL8 mCostGeomProhibited = 1e5;


cAVertReg::cAVertReg(eTypeSom aType,const  std::vector<cPt2di> & aVPts,int aNum) :
    mType       (aType),
    mVPts       (aVPts),
    mNumV       (aNum),
    mNumCC0     (-1),
    mNumMaxLoc  (-1),
    mProvMaxLoc  (-1),
    mXC          (-1.0),
    mBox         (cBox2di::FromVect(aVPts,true))  // true => can be empty
{
}


bool cAVertReg::IsBorderExt() const
{
    return mType == eTypeSom::eBordeExt;
}
bool cAVertReg::IsStd() const
{
    return mType == eTypeSom::eStd;
}

const cPt2di& cAVertReg::PtApp(int aK) const
{
   return mPath.at(mIndApP.at(aK));
}

tSeg2dr  cAVertReg::Tg1() const
{
   return tSeg2dr(ToR(PtApp(0)),ToR(PtApp(1)));
}

tSeg2dr  cAVertReg::Tg2() const
{
    size_t aNb = mIndApP.size();
   return tSeg2dr(ToR(PtApp(aNb-1)),ToR(PtApp(aNb-2)));
}
//t
//tSeg2dr  Tg2() const;


cAOrientReg::cAOrientReg() :
    mPt        (0,0),
    mCostGeom  (mCostGeomProhibited)
{
}

cASymReg::cASymReg() :
    mCostGeom  (mCostGeomProhibited),
    mMaxVal (-1e10),
    mInMST  (false)
{
}

void SetEdgeForce(tEdgeReg * anE)
{
    anE->AttrSym().mMaxVal= mMaxValForced;
    anE->AttrSym().mCostGeom = mCostGeomForced;

    anE->AttrOriented().mCostGeom = mCostGeomForced;
    anE->AttrInvOr().mCostGeom = mCostGeomForced;
}


static constexpr int LabelRegBorder= 0;
static constexpr int LabelRegUnInit = -1;
static constexpr int LabelRegMin= -2;






static constexpr tREAL8 mMarginTouchC = 3.0;
static constexpr tREAL8 mMarginTouchBordVert = 20.0;




class  cWeigtGrReg : public cAlgo_ParamVG<tGraphReg>
{
   public :
       cWeigtGrReg (bool isModeGeom) :
           mIsModeGeom (isModeGeom)
       {

       }
       tREAL8 WeightEdge(const    tEdgeReg & anEdge) const override
       {
           return   mIsModeGeom  ?
                    anEdge.AttrSym().mCostGeom :
                    (- anEdge.AttrSym().mMaxVal);
       }
       bool mIsModeGeom;


};

class  cGrForestNoBorder : public cAlgo_ParamVG<tGraphReg>
{
   public :

     cGrForestNoBorder (bool isModeGeom) :
        mIsModeGeom (isModeGeom)
     {

     }
     bool   InsideVertex(const  tVertex &aV) const override {return !aV.Attr().IsBorderExt() ;}

     bool   InsideEdge(const    tEdge & anE) const override
     {
         if  (! anE.AttrSym().mInMST)
             return false;
         if (mIsModeGeom)
             return anE.AttrSym().mCostGeom < mCostGeomProhibited;
         else
             return  true;
     }

     bool mIsModeGeom;

};


cPt2dr cAppliNewFrange::PointCorrectIntegral(const cPt2di& aPt) const
{
    return cPt2dr(aPt.x()-mDImIntegr->GetV(aPt.y()),aPt.y());
}

tREAL8 cAppliNewFrange::CostGeomP1P2(const cAVertReg *anAttr1, const cAVertReg *anAttr2)  const
{
    cPt2di aP1 = anAttr1->mP1;
    cPt2di aP2 = anAttr2->mP2;
    cPt2dr aV12 = ToR(aP2-aP1);

    tREAL8 aDistE = Norm2(aP1-aP2);

    // really too far,
    if (aDistE > mIntervMax)
        return mCostGeomProhibited;

    tREAL8 aX1Curv = aP1.x() -mDImIntegr->GetV(aP1.y());
    tREAL8 aX2Curv = aP2.x() -mDImIntegr->GetV(aP2.y());
    tREAL8 aDifXCurv = std::abs(aX1Curv-aX2Curv);

    if (aDifXCurv > mIntervMax)
        return mCostGeomProhibited;

    cSegment2DCompiled<tREAL8> aS1 (anAttr1->Tg1());
    cSegment2DCompiled<tREAL8>  aS2 (anAttr2->Tg2());

    cPt2dr aDir1 = - aS1.Tgt();
    cPt2dr aDir2 = aS2.Tgt();
    cPt2dr aDir12 = VUnit(aV12);

    tREAL8 aDifTeta = (Norm2(aDir1-aDir12) + Norm2(aDir2-aDir12) + Norm2(aDir1-aDir2)) / 3.0;

    if (aDifTeta>1)
        return mCostGeomProhibited;

   // tREAL8 aDLine1 = aS1.DistLine(ToR(aP2));
   // tREAL8 aDLine2 = aS2.DistLine(ToR(aP1));


    return      (0.2+aDifTeta)
            *   (
                     //aDLine1+aDLine2
                     aDistE
                   + aDifXCurv
                )
     ;
}





/*
tREAL8 cAppliNewFrange::CostGeomP1P2(const cAVertReg *anAttr1, const cAVertReg *anAttr2)  const
{
   // we connect the two closest point
   if (!OkGeomP1P2(anAttr1,anAttr2))
        return mCostGeomProhibited;

   tRE

   return mCostGeomProhibited;
}
*/


void cAppliNewFrange::ComputeArribute(cAVertReg& anAttr,cDataIm2D<tU_INT1> & aDIm)
{
    if ( ! anAttr.IsStd())
        return;

    // Compute extremities as diameters of Pts
     auto anExtres = Diameter(anAttr.mVPts,aDIm,Alloc8Neighbourhood(),1,2,&anAttr.mCenter,&anAttr.mPath);
     anAttr.mP1 = anExtres.first;
     anAttr.mP2 = anExtres.second;

     // Orientate
     {
         cPt2di aV12 =  anAttr.mP2- anAttr.mP1;

         // if rather vertical
         bool DoSwap = false;
         if (std::abs(aV12.y()) > std::abs(aV12.x()))
         {
                 DoSwap = (aV12.y()>0);
         }
         else
         {
            bool halfUp = (anAttr.mP1.y() +anAttr.mP2.y()) / 2.0 > mYC;
            DoSwap =  (halfUp == (aV12.x()>0));
         }
         // if swapped, path & approx must be coherent
         if (DoSwap)
         {
             std::swap(anAttr.mP1,anAttr.mP2);
             std::reverse(anAttr.mPath.begin(),anAttr.mPath.end());
         }
     }
     anAttr.mIndApP = PolynSimplify(anAttr.mPath,2.0,100,3.0,5.0);



    tREAL8 aXC = -1;
    cWhichMin<cPt2di,int> aWProjLine;

    // compute the intersection/porjection with midle line
    if (1)
    {
        // Compute average position of intersection with central axe
        cWeightAv<tREAL8,tREAL8> aAvgXc;

        for (const auto & aPix : anAttr.mVPts)
        {
            int aDiff = std::abs(aPix.y() - mYC);
            aWProjLine.Add(aPix,aDiff);
            if (aDiff == 0)
                aAvgXc.Add(1.0,aPix.x());
        }

        if (aAvgXc.Nb())
        {
           aXC = aAvgXc.Average();
        }
    }

    // try the most obvious case if it cross the midlle line, and is close enouh to a X-Frange, affect it
    if (aXC>=0)
    {
        anAttr.mXC =  aXC;

        for (size_t aK=0 ; aK<mVMaxLoc.size() ; aK++)
        {
            const cMaxLocSum & aMLS = mVMaxLoc.at(aK);
            if (
                (aMLS.mX0 < anAttr.mXC) && ( anAttr.mXC<aMLS.mX1)
                && (std::abs(anAttr.mXC-aMLS.mXC) < mIntervMax / 4.0)
                )
            {
                MMVII_INTERNAL_ASSERT_always(anAttr.mNumMaxLoc==-1,"Multiple init mNumMaxLoc ");
                anAttr.mNumMaxLoc = aK;
                anAttr.mProvMaxLoc  = 0;
            }
        }
    }

    // second tentative, try point close enough
    if (anAttr.mProvMaxLoc<0)
    {
        cPt2di aPPL = aWProjLine.IndexExtre();
        if (std::abs(mDImTgt->GetV(aPPL.y()))<1.0)
        {
            tREAL8 aXCorr = aPPL.x() - mDImIntegr->GetV(aPPL.y());
            for (size_t aK=0 ; aK<mVMaxLoc.size() ; aK++)
            {
                const cMaxLocSum & aMLS = mVMaxLoc.at(aK);
                if (
                       (std::abs(aXCorr-aMLS.mXC) < mIntervMax / 4.0)
                    )
                {
                    MMVII_INTERNAL_ASSERT_always(anAttr.mNumMaxLoc==-1,"Multiple init mNumMaxLoc ");
                    anAttr.mNumMaxLoc = aK;
                    anAttr.mProvMaxLoc  = 1;
                }
            }
        }
    }


  /*  if (anAttr.mProvMaxLoc  <0)
    {
        cWhichMin<cPt2di,int> aClosest;
        for (const auto & aML : mVMaxLoc)
        {
            for (const auto & aPix : anAttr.mVPts)
            {
                cPt2di aDif =  aPix-aML.mPt;
                aClosest.Add(aPix,SqN2(cPt2dr(aDif.x(),aDif.y()*aMulY)));
            }
        }
    }*/

    // Try to detect line that touch the left
    if (anAttr.mProvMaxLoc <0)
    {
        cWhichMin<cPt2di,int> aLeftest;
        for (const auto & aPix : anAttr.mVPts)
        {
            aLeftest.Add(aPix,aPix.x());
        }
        cPt2di aProj = aLeftest.IndexExtre();
        // int aDif = std::abs(aProj.y()-mYC);
        tREAL8 aXCorr = aProj.x() - mDImIntegr->GetV(aProj.y());

        if (    (aProj.x()<100.0)
            &&  (aXCorr < mVMaxLoc.at(0).mXC - mIntervMax/2.0)
            )

            //(aProj.x()<aDif) && (aProj.x()<100) && (aDif>50))
        {
            anAttr.mProvMaxLoc = 2;
            anAttr.mNumMaxLoc =  mVMaxLoc.size() ;
        }
    }

    // if we dont have intersection, try close enoug
 /*   if (anAttr.mProvMaxLoc  <0)
    {
         tREAL8 aThrs = Square(mIntervMax / 4.0);
         tREAL8  aMulY = 2.0;

         cWhichMin<cPt2di,int> aClosest;
         for (const auto & aML : mVMaxLoc)
         {
             for (const auto & aPix : anAttr.mVPts)
             {
                 cPt2di aDif =  aPix-aML.mPt;
                 aClosest.Add(aPix,SqN2(cPt2dr(aDif.x(),aDif.y()*aMulY)));
             }
         }

         if (aClosest.ValExtre()<aThrs)
         {
             aXC = aClosest.IndexExtre().x();
              anAttr.mProvMaxLoc  = 1;
         }
    }*/




    anAttr.mUpSide = anAttr.mBox.P1().y()     < mYC - mMarginTouchC;
    anAttr.mBottomSide = anAttr.mBox.P0().y() > mYC + mMarginTouchC;
    anAttr.mNoCross = anAttr.mUpSide || anAttr.mBottomSide;

    cPt2di aP1 = anAttr.mP1;
    cPt2di aP2 = anAttr.mP2;

    tREAL8 aXMin = std::min(aP1.x(),aP2.x());
    tREAL8 aXMax = std::max(aP1.x(),aP2.x());

    anAttr.mTouchLeft = aXMin < mMarginTouchBordVert ;
    anAttr.mTouchRight = aXMax > mSzRed.x() - mMarginTouchBordVert ;

}




void cAppliNewFrange::PostProcEdge(tEdgeReg& anE)
{
    int aN1 = anE.VertexInit().Attr().mNumMaxLoc;
    int aN2 = anE.Succ().Attr().mNumMaxLoc;


    if ((aN1>=0) && (aN2>=0))
    {
          if (aN1==aN2)
          {
               anE.AttrSym().mMaxVal = mMaxValForced;
          }
          else
          {
              anE.AttrSym().mMaxVal = mMaxValProhibited;
          }
    }
}

void cAppliNewFrange::PostProcVertex(tVertexReg& aV)
{
    if ( aV.Attr().IsBorderExt())
        return;

    if (aV.Attr().mNoCross)
    {
        bool aTL = aV.Attr().mTouchLeft;
        bool aTR = aV.Attr().mTouchRight;

        if (aTL)
        {
           for ( auto & anE : aV.EdgesSucc())
           {
               if (anE->Succ().Attr().IsBorderExt())
               {
                    anE->AttrSym().mMaxVal = mMaxValForced;
               }
               else
               {
                   if (aTL && aTR)
                   {
                        anE->AttrSym().mMaxVal = mMaxValProhibited;
                   }
               }
           }
        }
    }
}


void cAppliNewFrange::MakeImageRegion()
{
   // StdOut()  << " MakeImageRegionMakeImageRegionMakeImageRegion " << mGrReg.AllVertices().size()  << "\n";
    // compute the diameters
    {
        cIm2D<tU_INT1>  aImMarq(mSzRed,nullptr,eModeInitImage::eMIA_Null);
        for (const auto & aVert : mGrReg.AllVertices())
        {
            ComputeArribute(aVert->Attr(),aImMarq.DIm());

        }

      //  StdOut() << "SUM MARQ IM" <<  SumIm(aImMarq.DIm()) << "\n";


    }



    mImRegion =  cIm2D<tINT2> (mSzRed,nullptr,eModeInitImage::eMIA_Null);
    mDImRegion = & mImRegion.DIm();

    mDImRegion->InitCste(LabelRegUnInit);
    mDImRegion->InitBorder(LabelRegBorder);

    // Make barrier with min region
    for (const auto& aPix: *mDImRegion)
    {
        if (LabelIsOn(aPix,eTypePix::eMinTens))
            mDImRegion->SetV(aPix,LabelRegMin);
    }

    // init the vec point to propagate wiht max line
    std::vector<cPt3di> aVPCur;
    for (const auto & aVert : mGrReg.AllVertices())
    {
        const auto & aCC = aVert->Attr();

        for (const auto& aPt: aCC.mVPts)
        {
            mDImRegion->SetV(aPt,aCC.mNumV);
            aVPCur.push_back(TP3z(aPt,aCC.mNumV));
        }
    }
    const auto & aVP4 = Alloc4Neighbourhood();

    // do the propagation
    while (!aVPCur.empty())
    {
        std::vector<cPt3di> aNewVP;
        for (const auto & aP3 : aVPCur)
        {
            cPt2di aCurP2 = Proj(aP3);
            for (const auto & aN4 : aVP4)
            {
                cPt2di aNewP = aCurP2 + aN4;
                tINT2 & aV = mDImRegion->GetReference_V(aNewP);
                if ((aV==LabelRegUnInit) || (aV==LabelRegMin))
                {
                    if (aV==LabelRegUnInit)
                        aNewVP.push_back(TP3z(aNewP,aP3.z()));
                     aV = aP3.z();
                }
            }

        }
        aVPCur =    aNewVP;
    }



    //  construct the edge of the graph by identify pixels wih different label
    for (const auto & aPix1 : *mDImRegion)
    {
        for (const auto & aN4 : aVP4)
        {
            cPt2di aPix2 = aPix1 + aN4;
            if (mDImRegion->Inside(aPix2))
            {
                 int aNum1 = mDImRegion->GetV(aPix1);
                 int aNum2 = mDImRegion->GetV(aPix2);
                 if (      (aNum1!=aNum2) && (aNum1>=0) && (aNum2>=0)
                       // && (!LabelIsOn(aPix1,eTypePix::eMinTens))
                       // && (!LabelIsOn(aPix2,eTypePix::eMinTens))
                     )
                 {
                     if (
                                (aPix1.x() != mSzRed.x()-1) && (aPix2.x() != mSzRed.x()-1)
                             && (aPix1.y() != mSzRed.y()-1) && (aPix2.y() != mSzRed.y()-1)
                             && (aPix1.y() != 0           ) && (aPix2.y() != 0           )
                        )
                     {
                         tVertexReg & aS1 = mGrReg.VertexOfNum(aNum1);
                         tVertexReg & aS2 = mGrReg.VertexOfNum(aNum2);

                         tEdgeReg* anE = mGrReg.AddEdge(aS1,aS2,cAOrientReg(),cAOrientReg(),cASymReg(),true);

                         tREAL8 aCost = std::min(mDImRedBlur->GetV(aPix1),mDImRedBlur->GetV(aPix2));

                         if (aCost > anE->AttrSym().mMaxVal)
                         {
                             anE->AttrSym().mMaxVal = aCost;
                             anE->AttrOriented().mPt = aPix1;
                             anE->AttrInvOr().mPt = aPix2;
                         }
                     }
                 }
            }
        }
    }

    // Connect the Agrg region to master
    for (const auto& aV : mVAgrFrg)
    {
        tEdgeReg* anE = mGrReg.AddEdge(*aV,*mVirtExtVert,cAOrientReg(),cAOrientReg(),cASymReg());
        SetEdgeForce(anE);
//        anE->AttrSym().mMaxVal= mMaxValForced;
    }
    // Connect the verectex affected to a region to the agreg
    for (const auto& aV : mGrReg.AllVertices())
    {
         int aNumML = aV->Attr().mNumMaxLoc;
         if ((aV->Attr().IsStd()) && (aNumML>=0))
         {
             tEdgeReg* anE = mGrReg.AddEdge(*aV,*mVAgrFrg.at(aNumML),cAOrientReg(),cAOrientReg(),cASymReg());
             SetEdgeForce(anE);
//             anE->AttrSym().mMaxVal= mMaxValForced;
         }
    }


    bool isModeGeom = false;
    // Apply different heuristic rules to modiy cost
    /*for (auto & anE : mGrReg.AllEdges_DirInit())
    {
         PostProcEdge(*anE);
    }*/
    /*for (auto & aV : mGrReg.AllVertices())
    {
         PostProcVertex(*aV);
    }*/


    if (isModeGeom)
    {
       for (auto & anE : mGrReg.AllEdges_DirInit())
       {
             tVertexReg & aV1 = anE->VertexInit();
             tVertexReg & aV2 = anE->Succ();
             if (aV1.Attr().IsStd()&&aV2.Attr().IsStd())
             {
                tREAL8 aCost12 = CostGeomP1P2(&aV1.Attr(),&aV2.Attr());
                tREAL8 aCost21 = CostGeomP1P2(&aV2.Attr(),&aV1.Attr());
                anE->AttrOriented().mCostGeom = aCost12;
                anE->AttrInvOr().mCostGeom = aCost21;
                anE->AttrSym().mCostGeom = std::min(aCost12,aCost21);
             }

       }
    }
    //  compute the minimal forest corresponding to cost mMaxVal
    {
        cVG_Forest<tGraphReg>   aForest(mGrReg);
        cWeigtGrReg  aWGr(isModeGeom);
        cAlgoSP<tGraphReg> anAlgo;
        anAlgo.MinimumSpanningForest(aForest,mGrReg,mGrReg.AllVertices(),aWGr);
        for (const auto & aTree : aForest.VTrees())
        {
            for (const auto & anEdge : aTree.Edges())
                anEdge->AttrSym().mInMST = true;
        }
    }

    //  --- Compute connected component

    {
       cGrForestNoBorder aGrFNB(isModeGeom);
       cAlgoCC<tGraphReg>  anAlgo;
       std::list<std::vector<tVertexReg *>> aLCC = anAlgo.All_ConnectedComponent(mGrReg,aGrFNB);
       int aNumCC =0;

       for ( auto & aCC : aLCC)
       {
           for ( auto & aV : aCC)
           {
               aV->Attr().mNumCC0 = aNumCC;
              // StdOut() << "CC00=" << aV->Attr().mNumCC0 << "\n" ;
           }
          aNumCC++;
       }
    }
   // StdOut() << " VTREE " << aForest.VTrees().size() << "\n";
}

};  // NS_FrangesDetect



}; //  namespace MMVII
