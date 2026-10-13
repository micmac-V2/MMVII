//#include "MMVII_PCSens.h"

// Test commit

#include "cAppliNewFranges.h"




namespace MMVII
{

namespace NS_FrangesDetect
{

/* **************************************** */
/*           eTypeSom                       */
/* **************************************** */

eTypeSom TypeSomStd(const bool isHor,const bool isP1)
{
    if (isHor)
    {
        if (isP1)
            return eTypeSom::eH1;
        else
            return eTypeSom::eH2;
    }

    if (isP1)
        return eTypeSom::eV1;

    return eTypeSom::eV2;
}

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


}; //  mamespace NS_FrangesDetect
}; //  namespace MMVII

