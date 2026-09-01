#include "MMVII_Ptxd.h"
#include "MMVII_Image2D.h"
#include "MMVII_Geom2D.h"



namespace MMVII
{



/* =============================================== */
/*                                                 */
/*                  cQuadTreeCell                  */
/*                                                 */
/* =============================================== */
cQuadTreeCell::cQuadTreeCell(cPixBox<2> aArea, int aLevel):
    mArea(aArea), mLevel(aLevel)
{
    mSubs.reserve(4);
}

cQuadTreeCell::cQuadTreeCell(cPixBox<2> aAreaInit):
mArea(aAreaInit), mLevel(1)
{
    auto aULPoint = mArea.P0();
    auto aLRPoint = mArea.P1();
    auto aCCPoint = (aULPoint+aLRPoint)/2;
    auto aUCPoint = cPt2di(aCCPoint.x(), mArea.P0().y());
    auto aLCPoint = cPt2di(aCCPoint.x(), mArea.P1().y());
    mSubs.push_back(cQuadTreeCell(cPixBox<2>(aULPoint,aLCPoint),mLevel+1));
    mSubs.push_back(cQuadTreeCell(cPixBox<2>(aUCPoint,aLRPoint),mLevel+1));
}


bool cQuadTreeCell::IsLeaf() const
{
    return mSubs.empty();
}

const cPixBox<2> & cQuadTreeCell::GetArea() const
{
    return mArea;
}

int cQuadTreeCell::Divide4(int aMinCellSz)
{
    if (mSubs.empty())
    {
        auto aULPoint = mArea.P0();
        auto aLRPoint = mArea.P1();
        //std::cout<<"  Split "<<aULPoint<<" "<<aLRPoint;
        if ((aLRPoint.x()-aULPoint.x()<aMinCellSz) || (aLRPoint.y()-aULPoint.y()<aMinCellSz))
        {
            //std::cout<<" too small\n";
            return 0; // do not split too small
        }
        auto aCCPoint = (aULPoint+aLRPoint)/2;
        auto aUCPoint = cPt2di(aCCPoint.x(), mArea.P0().y());
        auto aCLPoint = cPt2di(mArea.P0().x(), aCCPoint.y());
        auto aCRPoint = cPt2di(mArea.P1().x(), aCCPoint.y());
        auto aLCPoint = cPt2di(aCCPoint.x(), mArea.P1().y());
        mSubs.push_back(cQuadTreeCell(cPixBox<2>(aULPoint,aCCPoint),mLevel+1));
        mSubs.push_back(cQuadTreeCell(cPixBox<2>(aUCPoint,aCRPoint),mLevel+1));
        mSubs.push_back(cQuadTreeCell(cPixBox<2>(aCLPoint,aLCPoint),mLevel+1));
        mSubs.push_back(cQuadTreeCell(cPixBox<2>(aCCPoint,aLRPoint),mLevel+1));
        //std::cout<<" ok\n";
        return mSubs.size();
    } else {
        // do nothing
        return 0;
    }
}



/* =============================================== */
/*                                                 */
/*                    cQuadTree                    */
/*                                                 */
/* =============================================== */
cQuadTree::cQuadTree(cDataIm2D<tREAL4> *aDepthIm, int aMinCellSz):
    mDepthIm(aDepthIm), mCurNbCell(1),
    mRootCell(*mDepthIm), mMinCellSz(aMinCellSz)
{

}

int cQuadTree::GetCurNbCell() const
{
    return mCurNbCell;
}

void cQuadTree::Split(int aTargetNbCell)
{
    std::cout<<"Split into "<<aTargetNbCell<<" cells...\n";
    mVLeafs.clear();
    float aLimit =  mDepthIm->MoyVal() / std::log2(aTargetNbCell) * 2.; // > to what we expect, to make room for variation
    while ((mCurNbCell<aTargetNbCell)&&(aLimit>1e-3))
    {
        //std::cout<<"Split step\n";
        // loop over all subcells
        mCurNbCell+=RecursiveDiv(&mRootCell, aLimit);
        aLimit *= 0.98; // change limit to add cells if needed, slowly to get close to asked number of cells
        //std::cout<<"CurNbCell: "<<mCurNbCell<<"\n";
    }

    RecursiveFillVLeafs(&mRootCell);
}

int cQuadTree::RecursiveDiv(cQuadTreeCell *aCell, float aLimit)
{
    int aNbCreated = 0;
    if (aCell->mSubs.empty())
    {
        // test if dist_max/level too big, sqrt on dist to limit a little the contrast between close and far?
        float aMax = mDepthIm->MaxVal(aCell->mArea);
        float aMin = mDepthIm->MinValNotNull(aCell->mArea);
        float aVal = sqrt(aMax * sqrt(1+sqrt((aMax-aMin)/(aMax+1.)))) / Square(aCell->mLevel);
        //std::cout<<"  aVal="<<aVal<<" aLimit="<<aLimit<<"\n";
        if (aVal > aLimit)
            aNbCreated += aCell->Divide4(mMinCellSz) - 1; // lost 1 leaf, made n new leaves
    }
    for (auto & aSubCell: aCell->mSubs)
        aNbCreated += RecursiveDiv(&aSubCell, aLimit);
    return aNbCreated;
}

void cQuadTree::RecursiveFillVLeafs(cQuadTreeCell * aCell)
{
    if (aCell->IsLeaf())
        mVLeafs.push_back(aCell);
    else
    {
        for (auto & aSubCell: aCell->mSubs)
            RecursiveFillVLeafs(&aSubCell);
    }
}

const std::vector<const cQuadTreeCell*> cQuadTree::GetVLeafs() const
{
    return mVLeafs;
}

};
