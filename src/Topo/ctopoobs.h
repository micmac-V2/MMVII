#ifndef CTOPOOBS_H
#define CTOPOOBS_H

#include "MMVII_AllClassDeclare.h"
#include "MMVII_enums.h"
#include "MMVII_SysSurR.h"
#include "MMVII_Geom3D.h"

namespace MMVII
{
class cTopoObsSet;
class cTopoPoint;
class cBA_Topo;
class cTopoObs;

template <class Type> class cResidualWeighter;

/**
 * @brief The cTopoSigma class represents the precision of an observation.
 * it is computed from obs length:
 *  SigmaTotal = SigmaAbs + SigmaRel * length for distances
 *  SigmaTotal = SigmaAbs + SigmaRel / length for angles
 *  w = 1/s**2
 */
class cTopoSigma
{
public:
    cTopoSigma(); // for serialization
    cTopoSigma(tREAL8 aSigmaAbs, tREAL8 aSigmaRel);
    tREAL8 getTotalSigma(const cTopoObs &aTopoObs, tREAL8 aLength) const;
    tREAL8 getTotalWeight(const cTopoObs &aTopoObs, tREAL8 aLength) const;
    tREAL8 mSigmaAbs;
    tREAL8 mSigmaRel;
};
void AddData(const cAuxAr2007 & anAuxInit, cTopoSigma & aTopoSigma);


/**
 * @brief The cTopoWeighter class records an obs cTopoSigmas
 * and get the total weight for least squares using current length
 */
class cTopoWeighter : public cResidualWeighter<tREAL8>
{
public :
    typedef std::vector<tREAL8>     tStdVect;

    cTopoWeighter(cTopoObs* aTopoObs);
    std::vector<tREAL8> WeightOfResidual(const tStdVect &) const override;
private :
    cTopoObs* mTopoObs;
};



/**
 * @brief The cTopoObs class represents an observation between several points.
 * It exists only in a cTopoObsSet because the obs may share parameters with other obs.
 */
class cTopoObs : public cMemCheck
{
    friend class cTopoObsSet;
    friend class cTopoData;
public:
    //~cTopoObs() { std::cout<<"delete topo obs "<<toString()<<std::endl; }
    std::string toString() const;
    eTopoObsType getType() const {return mType;}
    std::vector<int> getIndices() const;
    std::vector<tREAL8> getVals() const; //< for least squares (with rotation matrix if needed)
    const std::vector<tREAL8> & getMeasures() const { return mMeasures;} //< original measures
    std::vector<tREAL8> & getResiduals() { return mLastResiduals;} //< last residuals
    const std::vector<tREAL8> & getResiduals() const { return mLastResiduals;} //< last residuals
    cTopoWeighter& getWeights();
    const std::string & getPointName(size_t i) const { return mPtsNames.at(i); }
    const std::vector<std::string> & getPointNames() const { return mPtsNames; }
    const std::vector<cTopoSigma> & getTopoSigmas() const { return mTopoSigmas; }
    tREAL8 getLength() const;
    //std::vector<tREAL8> getResiduals(const cTopoComp *comp) const;
protected:
    cTopoObs(cTopoObsSet* set, cBA_Topo * aBA_Topo, eTopoObsType type, const std::vector<std::string> & ptsNames,
             const std::vector<tREAL8> & measures,  const std::vector<cTopoSigma> &aTopoSigmas);
    cTopoObs(const cTopoObs &) = delete;
    cTopoObs& operator=(const cTopoObs &) = delete;
    cTopoObsSet* mSet;//the set containing the shared parameters
    cBA_Topo * mBA_Topo;
    eTopoObsType mType;
    //std::vector<cTopoPoint*> mPts;
    std::vector<std::string> mPtsNames;
    std::vector<tREAL8> mMeasures;
    std::vector<cTopoSigma> mTopoSigmas;
    cTopoWeighter mWeights;
    std::vector<tREAL8> mLastResiduals;
};

};
#endif // CTOPOOBS_H
