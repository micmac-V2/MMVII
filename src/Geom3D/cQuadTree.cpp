#include "MMVII_Ptxd.h"
#include "MMVII_Image2D.h"
#include "MMVII_QuadTree.h"



namespace MMVII
{



/* =============================================== */
/*                                                 */
/*                  cQuadTreeCell                  */
/*                                                 */
/* =============================================== */
cQuadTreeCell::cQuadTreeCell(cPixBox<2> aArea, int aLevel):
    mArea(aArea), mLevel(aLevel), mValsComputed(false), mValMin(NAN), mValMax(NAN)
{
    mSubs.reserve(4);
}

cQuadTreeCell::cQuadTreeCell(cPixBox<2> aAreaInit):
mArea(aAreaInit), mLevel(1), mValsComputed(false), mValMin(NAN), mValMax(NAN)
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

int cQuadTreeCell::DivideHV(int aNbX, int aNbY)
{
    if (mSubs.empty())
    {
        auto aStep = cPt2di(mArea.P1() - mArea.P0())/aNbX;

        for (int aY=0;aY<aNbY;aY++)
            for (int aX=0;aX<aNbX;aX++)
            {
                auto aUL = cPt2di(mArea.P0().x() + aStep.x()*aX, mArea.P0().y() + aStep.y()*aY);
                auto aLR = cPt2di(mArea.P0().x() + aStep.x()*(aX+1), mArea.P0().y() + aStep.y()*(aY+1));
                mSubs.push_back(cQuadTreeCell(cPixBox<2>(aUL,aLR),mLevel+1));
            }
        return mSubs.size();
    } else {
        // do nothing
        return 0;
    }
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
    int aIter = 0;
    int aPrevNbCell = 0;

    // first split is in a number of cells depending on image ratio
    int aStartNumberCell = std::min(100,aTargetNbCell);
    int aStartNX = sqrt(aStartNumberCell * mDepthIm->SzX() / mDepthIm->SzY());
    int aStartNY = sqrt(aStartNumberCell * mDepthIm->SzY() / mDepthIm->SzX());
    mRootCell.DivideHV(aStartNX, aStartNY);

    while ((mCurNbCell<aTargetNbCell)&&(aIter<1000))
    {
        aPrevNbCell = mCurNbCell;
        //std::cout<<"Split step\n";
        // loop over all subcells
        mCurNbCell+=RecursiveDiv(&mRootCell, aLimit);
        aLimit *= 0.98; // change limit to add cells if needed, slowly to get close to asked number of cells
        if (mCurNbCell==aPrevNbCell)
            aLimit *= 0.8; // accelerate when not changing
        //std::cout<<"CurNbCell: "<<mCurNbCell<<"\n";
        aIter++;
    }
    std::cout<<"Got "<<mCurNbCell<<" cells in "<<aIter<<" iterations\n";

    RecursiveFillVLeafs(&mRootCell);
}

int cQuadTree::RecursiveDiv(cQuadTreeCell *aCell, float aLimit)
{
    int aNbCreated = 0;
    if (aCell->mSubs.empty())
    {
        // test if dist_max/level too big, sqrt on dist to limit a little the contrast between close and far?
        if (!aCell->mValsComputed)
        {
            aCell->mValMax = mDepthIm->MaxVal(aCell->mArea);
            aCell->mValMin = mDepthIm->MinValNotNull(aCell->mArea);
            aCell->mValsComputed = true;
        }
        float aVal = (sqrt(sqrt(aCell->mValMax)) * (1+(aCell->mValMax-aCell->mValMin)/(aCell->mValMax+.01))) / Square(Square(aCell->mLevel));
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
