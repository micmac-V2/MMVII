#ifndef _MMVII_SelecTripletsCoverage_H_
#define _MMVII_SelecTripletsCoverage_H_

#include "MMVII_Ptxd.h"

namespace MMVII
{

class cParamSelTriCov;        // parameters of the coverage/ B/H triplet selection
class cGridCoverTri;          // grid on the bounding box of the pair tie points, in image A
class cCandTriCov;            // one candidate third image, as given by the caller
class cSelecTripletsCoverage; // greedy selection of third images for one oriented pair

/// Parameters of the selection (defaults = MMV1 Martini)
class cParamSelTriCov
{
public :
    cParamSelTriCov();

    bool    mGreedy;     ///< if false, pass-through : every candidate is selected
    size_t  mNbMax;      ///< max number of third images selected for one pair
    tREAL8  mKappa;      ///< cap of the update in cells not covered by the selected image (MMV1 : 0.3)
    tREAL8  mGainThr;    ///< stop when the best normalised gain is below
    int     mNbCells;    ///< max number of cells per side of the grid
    tREAL8  mBsHLim;     ///< half-saturation of the B/H gain
    tREAL8  mAttenDens;  ///< saturation of the coverage term
};

/**  Square grid on the bounding box of the A-B tie points, in calibrated coordinates of image A.
 *   Cell side = max(width,height)/NbCells; points outside the box are clamped in border cells.
 */
class cGridCoverTri
{
public :
    cGridCoverTri(const std::vector<cPt2dr> & aVPts,int aNbCellsMax);

    size_t NbCells() const;                  ///< total number of cells
    size_t IndexOf(const cPt2dr &) const;    ///< index of the cell containing the point (clamped)
    std::vector<tREAL8> Count(const std::vector<cPt2dr> &) const; ///< number of points per cell

private :
    cPt2dr  mP0;     ///< lower corner of the box
    tREAL8  mStep;   ///< side of a cell
    int     mNbX;    ///< number of cells in x
    int     mNbY;    ///< number of cells in y
};

/// One candidate third image C for the pair (A,B)
class cCandTriCov
{
public :
    int                  mId;      ///< identifier given by the caller (e.g. num of image C)
    cPt3dr               mCenter;  ///< centre of C in the frame of A, scale |C_B|=1
    std::vector<cPt2dr>  mVPts3;   ///< triple tie points A-B-C, coordinates in image A
};

/**  Greedy selection of third images for one oriented pair (A,B), port of MMV1 NO_GenTripl.
 *   Gain of C = mean over cells of  D_C[k] * gamma_C * P[k], where D_C is the coverage of the
 *   triple points, gamma_C the B/H gain and P the density of the pair points. After each
 *   selection of S, the other gains are capped by D_C[k] * (gamma(C,S) D_S[k] + Kappa (1-D_S[k])).
 */
class cSelecTripletsCoverage
{
public :
    cSelecTripletsCoverage
        (
            const cParamSelTriCov &     aParam,
            const std::vector<cPt2dr> & aVPts2,    ///< A-B tie points, coordinates in image A
            const cPt3dr &              aCenterB,  ///< centre of B in the frame of A
            tREAL8                      aDepthH    ///< depth of the scene, same scale
            );

    void AddCandidate(const cCandTriCov &);

    /// Run the selection (once : modifies the gains) ; return (Id,Gain at selection) in order of selection
    std::vector<std::pair<int,tREAL8>> Select();


    /// Coverage term of a density d in [0,1] :  (T+1)/T * T d/(T d+1)
    static tREAL8 CoverOfDensity(tREAL8 aDens,tREAL8 aAttenDens);
    /// B/H gain : bsh/(bsh+Lim)
    static tREAL8 GainOfBsH(tREAL8 aBsH,tREAL8 aBsHLim);

private :
    /// internal state of one candidate
    class cCandInternal
    {
    public :
        int                  mId;
        cPt3dr               mCenter;  ///< centre of C in the frame of A (|C_B|=1) : lambda * t_AC, lambda from triple points
        std::vector<tREAL8>  mCover;  ///< D_C[k]
        std::vector<tREAL8>  mGain;   ///< G_C[k], decreases with selections
        bool                 mSelected;
    };

    tREAL8 GainBsH(const cPt3dr & aC1,const cPt3dr & aC2) const; ///< gamma between 2 centres
    tREAL8 TotalGain(const cCandInternal &) const;               ///< (1/NbCells) Sum G[k] P[k]
    void   UpdateGain(cCandInternal & aCand,const cCandInternal & aSel); ///< MMV1 UpdateCost

    cParamSelTriCov             mParam;
    cGridCoverTri               mGrid;
    std::vector<tREAL8>         mNb2;      ///< number of pair points per cell
    std::vector<tREAL8>         mDensP;    ///< P[k] = sqrt(Nb2/maxNb2)
    cPt3dr                      mCenterB;
    tREAL8                      mDepthH;
    std::vector<cCandInternal>  mVCand;
};


} // MMVII

#endif //_MMVII_SelecTripletsCoverage_H_
