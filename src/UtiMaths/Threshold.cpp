#include "MMVII_util.h"
#include "MMVII_Stringifier.h"

#include "MMVII_BundleAdj.h"
namespace MMVII
{

cVarThreshold::cVarThreshold(double aThresholdInit, double aThresholdFinal,
              float aPartShrinking, bool aIsLinear) :
    mThresholdInit(aThresholdInit), mThresholdFinal(aThresholdFinal),
    mPartShrinking(aPartShrinking), mIsLinear(aIsLinear)
{
    MMVII_INTERNAL_ASSERT_User((mPartShrinking>=0.)&&(mPartShrinking<=1.),
                               eTyUEr::eUnClassedError,"PartShrinking must be in [0,1], got "+ToStr(mPartShrinking));

    //mThresholdInit = (aThresholdInit<0) ? INFINITY : aThresholdInit;
    //mThresholdFinal = (aThresholdFinal<0) ? INFINITY : aThresholdFinal;
}

double cVarThreshold::GetThreshold(int aIter, int aNbMaxIter) const
{
    // spend last xx% of iter on final threshold
    int aNbIterFinalTh = aNbMaxIter*mPartShrinking;

    double aThreshold = mThresholdFinal;
    if (aNbIterFinalTh<2)
        return aThreshold;

    if (mIsLinear)
        aThreshold = mThresholdInit + (mThresholdFinal - mThresholdInit)*float(aIter)/(aNbIterFinalTh-1);
    else
    {
        if (aIter <= aNbIterFinalTh)
            //mThreshold = mThresholdInit + (mThresholdFinal - mThresholdInit)*float(aBA.Iter())/(aNbIterFinalTh);
            aThreshold = mThresholdFinal + (mThresholdInit - mThresholdFinal) * (1. + cos(M_PI * float(aIter)/aNbIterFinalTh ))/2. ;
    }

    std::cout << "up weighters, th="<<aThreshold<<"\n";
    if (aThreshold>10000)
        aThreshold = -1;
    return aThreshold;
}

bool cVarThreshold::HasChangedSincePreviousIter(int aIter, int aNbMaxIter) const
{
    if (mThresholdInit == mThresholdFinal)
        return false;
    int aNbIterFinalTh = aNbMaxIter*mPartShrinking;
    return (aIter <= aNbIterFinalTh);
}

};
