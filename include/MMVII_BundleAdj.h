#ifndef  _MMVII_BUNDLE_ADJ_H_
#define  _MMVII_BUNDLE_ADJ_H_

#include "MMVII_Sensor.h"

using namespace NS_SymbolicDerivative;

namespace MMVII
{
/**  Class for otimizing a model of camera  using 3d-2d correspondance and bundle adjustment.  Typically these
 *   corresponance can be synthetic ones coming from another camera or "true" one comming from GCP.
 *  . It can be used in, two scenario :
 *
 *    -(1) primary test/bench  on functionnality to do BA
 *    -(2)
 *        (2.a)   conversion between calibration (format/model ...)
 *        (2.b)   comparison of calibrations (to come)
 *    -(3)  Optimization after space resection  , 11 parameters
 *    -(4)  Optimization after space resection
 *
 *    In first case we create artifcially difficult conditions (randomize the initial pose, let free the perspective center).
 *
 *    In the second case,  we use as much information we have : init with identity, and froze the position center
 *
 *    Also the general bundle class could be used, and would be preferable for avoiding code duplication, we maintain it as:
 *
 *       - it's a tiny example that can be used as a tuto
 *       - the case is common, and it easier to use and probably faster that general bundle class
 *       - as it's tiny there is not so much code duplication
 *
 */

class cCorresp32_BA
{
       public :
         typedef cIsometry3D<tREAL8>   tPose;

         cCorresp32_BA
         (
              cSensorImage       *,
              const cSet2D3D &
         );
         ~cCorresp32_BA();
         void OneIteration();
         const cSet2D3D  & SetCorresp() const {return   mSetCorresp;}
         cSet2D3D  & SetCorresp()  {return   mSetCorresp;}

         void SetFrozenVarOfPattern(const std::string & aPattern);

         cResolSysNonLinear<double> &   Sys();   ///< Solver
       protected :

         cSensorImage*                      mSensor;
         // When using conversion for real application, these two variable will be set to true, because we have no interest to hide
         // information. BTW in bench mode, we put the system in more difficult condition, to check that we all the same
         // get to the good solution (but a litlle slower)
         bool                               mFGC; // HardConstrOnGCP if true the 3d point are frozen
         bool                               mCFix; // Center Fix : if true center of rotation is frozen

         cSet2D3D                           mSetCorresp;  ///<  Set of 2D-3D correspondance
         int                                mSzBuf;   ///<  Sz Buf for calculator
         cCalculator<double> *              mEqColinearity;  ///< Colinearity equation
         cSetInterUK_MultipeObj<double>     mSetInterv;   ///< coordinator for autom numbering
         cResolSysNonLinear<double> *       mSys;   ///< Solver
};

/**
 * @brief The cVarThreshold class computes a threshold interpolated between
 * aThresholdInit and aThresholdFinal, depending on the current iteration.
 * The standard interpolation is in cosinus (slowly changing at start and end, quickly in the middle).
 * The interpolation can happen only on a part of the iterations:
 * aPartShrinking = 0.75 means that the last 25% of iterations are all at aThresholdFinal.
 * if isLinear, the interpolation is linear instead of cos.
 *
 * a value >10000 is transformed into -1 (meaning no limit for weighters)
 */
class cVarThreshold
{
public:
    cVarThreshold(double aThresholdInit, double aThresholdFinal,
                  float aPartShrinking = 0.6, bool aIsLinear=false);
    double GetThreshold(int aIter, int aNbMaxIter) const;
    bool HasChangedSincePreviousIter(int aIter, int aNbMaxIter) const;
protected:
    double mThresholdInit, mThresholdFinal;
    float mPartShrinking;
    bool mIsLinear;
};

};  //  MMVII

#endif  //  _MMVII_BUNDLE_ADJ_H_
