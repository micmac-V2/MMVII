//#include "MMVII_PCSens.h"

// Test commit

#include "cAppliNewFranges.h"




namespace MMVII
{

namespace NS_FrangesDetect
{


class cTensAVertex
{
    public :
        cPt2di   mPix;
        tREAL8   mCost;
};


struct cTensASym
{
   public :
      cTensASym(tREAL8 aCost);

      tREAL8 mCost;
};

struct cTensAOriented
{
    public :
       cTensAOriented(const tREAL8 aCost) ;
       tREAL8 mCost;
};


typedef cVG_Vertex<cTensAVertex,cTensAOriented,cTensASym> tTensVertex;
typedef cVG_Edge<cTensAVertex,cTensAOriented,cTensASym>   tTensEdge;
typedef cVG_Graph<cTensAVertex,cTensAOriented,cTensASym>  tTensGraph;




//typedef cVG_Edge<cAVVertex,cAVEdgeOriented,cAVEdgeSym>  tEdgeVFrange;
//typedef cVG_Graph<cAVVertex,cAVEdgeOriented,cAVEdgeSym> tGraphVFrange;



class cCOlTensVersex : public cMemCheck
{
    public :
      int mX0;
      int mX1;
      int mXC;
      int mY;

      cPt2di P0() const {return cPt2di(mX0,mY);}
      cPt2di P1() const {return cPt2di(mX1,mY);}
      cPt2di PC() const {return cPt2di(mXC,mY);}

};


//  int aDy = 1; // -1
//  int aYLim = mSzRed.y() ; // -1

void cAppliNewFrange::MakeGrapheTensor(int aXC0,int aStepX,int aStepY, cRGBImage* aPtrIm)
{
    tTensGraph * aNewGr = new tTensGraph();
    std::vector<cCOlTensVersex> aVCol;


    cCOlTensVersex aCol;
    for (int anY=mYC ; (anY>=0) && (anY<mSzRed.y()) ; anY+= aStepY )
    {
        int aDeltaY = std::abs(anY-mYC);
       // int aDXCurve = mDImIntegr->GetV(anY);
        int aXC = aXC0 + mDImIntegr->GetV(anY);;
        int aTgt = std::abs(mDImTgt->GetV(anY));
        aCol.mY = anY;
        int anIncX0 = std::min(round_ni(mIntervMax/2),std::max(round_ni(mIntervMax/4),aDeltaY+aTgt*10));
        int anIncX1 = std::min(round_ni(mIntervMax*1.5),std::max(round_ni(mIntervMax/4),aDeltaY*3+aTgt*20));

        aCol.mX0 = round_to(round_ni(aXC-anIncX0),aStepX);
        aCol.mX1 = round_to(round_ni(aXC+anIncX1),aStepX);

        while (aCol.mX0<0) aCol.mX0 += aStepX;
        while (aCol.mX1>=mSzRed.x()) aCol.mX1 -= aStepX;

        aCol.mXC = round_to(aXC ,aStepX);

        if (aCol.mX0<aCol.mX1)
            aVCol.push_back(aCol);
    }
    if (aPtrIm)
    {
        for (size_t aK=1 ;aK<aVCol.size(); aK++ )
        {
            auto aCol0 = aVCol.at(aK-1);
            auto aCol1 = aVCol.at(aK);

            aPtrIm->DrawLine(ToR(aCol0.P0()),ToR(aCol1.P0()),cRGBImage::Blue,1.0);
            aPtrIm->DrawLine(ToR(aCol0.P1()),ToR(aCol1.P1()),cRGBImage::Red,1.0);
            aPtrIm->DrawLine(ToR(aCol0.PC()),ToR(aCol1.PC()),cRGBImage::Green,1.0);

        }
        /*
        aPtrIm->SetRGBPix(cPt2di(aXC+aDXCurve,anY),cRGBImage::Green);

        aPtrIm->SetRGBPix(cPt2di(aX0,anY),cRGBImage::Blue);
        aPtrIm->SetRGBPix(cPt2di(aX1,anY),cRGBImage::Red);*/

    }
    delete aNewGr;

}

void cAppliNewFrange::MakeGrapheTensor(const cIm2D<tU_INT1> anIm)
{
    cRGBImage aRGB(anIm.DIm());

    int aXC = mVMaxLoc.at(0).mXC;
    MakeGrapheTensor(aXC,1,1,&aRGB);
    MakeGrapheTensor(aXC,1,-1,&aRGB);

    GenVisu(aRGB,"GrTens");
}




#if (0)
/* **************************************** */
/*           cASom                          */
/* **************************************** */

cAVVertex::cAVVertex(const cPt2di & aPt,const cPt2dr & aTgt,eTypeSom aType) :
    mPt       (aPt),
    mTgt      (aTgt),
    mType     (aType)
{
}

/* **************************************** */
/*           cASom                          */
/* **************************************** */


cAVEdgeSym::cAVEdgeSym(cConnComp* aCC):
    mCC (aCC)
{
}

/* **************************************** */
/*           cAEdgeOriented                 */
/* **************************************** */

cAVEdgeOriented::cAVEdgeOriented(tREAL8 aCost) :
    mCost (aCost)
{
}

/* **************************************** */
/*           cAppliNewFrange                */
/* **************************************** */

void cAppliNewFrange::GraphCCHoriz()
{

}


void cAppliNewFrange::GraphCreate()
{
    for (size_t aKC=0 ; aKC<mVectCC.size() ; aKC++)
    {
        cConnComp* aCC = mVectCC.at(aKC);

        if (aCC->mIsOk)
        {
            const tSeg2dr & aSeg = aCC->mSeg;

            cAVVertex aA1(ToI(aSeg.P1()),aCC->mTg1,TypeSomStd(aCC->mIsHor,true));
            cAVVertex aA2(ToI(aSeg.P2()),aCC->mTg2,TypeSomStd(aCC->mIsHor,false));
            cAVEdgeSym anAS(aCC);
            cAVEdgeOriented anA12(0.0);
            cAVEdgeOriented anA21(0.0);

            tVertexVFrange* aV1 = mGrVect.NewVertex(aA1);
            tVertexVFrange* aV2 = mGrVect.NewVertex(aA2);
            mGrVect.AddEdge(*aV1,*aV2,anA12,anA21,anAS);
        }
    }

    StdOut() << " NbVertex="  << mGrVect.NbVertex() << "\n";
}
//cASom::cASom(const cPt2di & aPt,const cPt2dr & aTgt,eTypeSom aType) :


void cAppliNewFrange::GraphCompute()
{
    GraphCreate();

    GraphCCHoriz();
}
#endif

}; //  mamespace NS_FrangesDetect
}; //  namespace MMVII

