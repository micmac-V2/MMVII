#ifndef  _MMVII_QUADTREE_H_
#define  _MMVII_QUADTREE_H_

#include "MMVII_Image2D.h"

namespace MMVII
{



class cQuadTreeCell
{
    friend class cQuadTree;
public:
    bool IsLeaf() const;
    const cPixBox<2> & GetArea() const;
protected:
    cQuadTreeCell(cPixBox<2> aArea, int aLevel);
    cQuadTreeCell(cPixBox<2> aAreaInit); //for inital level: splits area to get good proportions
    int Divide4(int aMinCellSz); // return number of cells created
    int DivideHV(int aNbX, int aNbY); // return number of cells created
    cPixBox<2> mArea;
    int mLevel;
    std::vector<cQuadTreeCell> mSubs;

    bool mValsComputed;
    float mValMin; // info saved to save time...
    float mValMax;
};



class cQuadTree
{
public:
    cQuadTree(cDataIm2D<tREAL4> * aDepthIm, int aMinCellSz=4);
    void Split(int aTargetNbCell);
    int GetCurNbCell() const;
    const std::vector<const cQuadTreeCell*> GetVLeafs() const;

protected:
    cDataIm2D<tREAL4> * mDepthIm;
    int mCurNbCell;
    cQuadTreeCell mRootCell;
    int RecursiveDiv(cQuadTreeCell * aCell, float aLimit); // returns number of created cells
    void RecursiveFillVLeafs(cQuadTreeCell * aCell);
    std::vector<const cQuadTreeCell*> mVLeafs;
    int mMinCellSz;
};



};

#endif  //  _MMVII_QUADTREE_H_
