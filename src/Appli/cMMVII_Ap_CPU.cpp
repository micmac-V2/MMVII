#include "cMMVII_Appli.h"
#include "MMVII_Sys.h"
#include <atomic>
#include <unordered_map>
#include <algorithm>
#include <vector>

namespace MMVII
{

std::string timePointAsString(const std::chrono::system_clock::time_point& tp)
{
    std::time_t t = std::chrono::system_clock::to_time_t(tp);
    std::string ts = std::ctime(&t);
    ts.resize(ts.size()-1);
    return ts;
}

// 3600.0 * 24.0

/***********************************/
/*                                 */
/*           cMMVII_Ap_CPU         */
/*                                 */
/***********************************/

cMMVII_Ap_CPU::cMMVII_Ap_CPU() :
   mT0             (std::chrono::system_clock::now())     ,
   mPid            (mmvii_GetPId())   ,
   mNbProcSystem   (mmvii_NbProcSys()),
   mNbProcAllowed  (mNbProcSystem),
   mMulNbInMk      (8.0),
   mTimeSegm       (this)
{
  // Very tricky and dirty, but I dont have courage for now to understand time/clok in C++
  // To change one day ....
  {
     const auto p0 = std::chrono::time_point<std::chrono::system_clock>{};  // Epoch 1/1/70
     std::chrono::duration<double, std::milli> fp_ms = mT0 -p0;  // duration from epoch
     double aT = fp_ms.count()/1000.0;  // In millisecond
     double  aSecPerDay = 24 * 3600;   // 24h/Day, 3600 sec/h ....
     aT = aT - (365.24219 *49 )  * aSecPerDay;  // Make epoch more or less to 1/1/2019 , use gregorian year in day
     aT = aT / aSecPerDay;  ///
     aT += 0.76736 -  0.816726 ;  // Experimental difference !!
     int aNbDay = round_down(aT);
     aT = aT - aNbDay;  // take fractionnal part
     aT = aT * aSecPerDay ;  //
     int aSec = round_down(aT);
     int aMili10 = std::min(9999,round_ni(1e4*(aT-aSec)));
     mStrIdTime =   ToStr(aNbDay,4)  + "_" + ToStr(aSec,5) +"_" + ToStr(aMili10,4) ;
  }
}

std::string    cMMVII_Ap_CPU::StrDateBegin() const    {return  timePointAsString(mT0);}
std::string    cMMVII_Ap_CPU::StrDateCur() const    {return  timePointAsString(std::chrono::system_clock::now());}
const std::string  &  cMMVII_Ap_CPU::StrIdTime() const  {return mStrIdTime;}


cTimerSegm&  cMMVII_Ap_CPU::TimeSegm() {return mTimeSegm;}

double cMMVII_Ap_CPU::SecFromT0() const
{
    tTime aT1 = std::chrono::system_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = aT1-mT0;
    return fp_ms.count()/1000.0;

}

/***********************************/
/*                                 */
/*          cTimeSequencer         */
/*                                 */
/***********************************/

cTimeSequencer::cTimeSequencer(double aPeriod):
    mPeriod   (aPeriod),
    mLastime  (cMMVII_Appli::CurrentAppli().SecFromT0())
{
}

bool cTimeSequencer::ItsTime2Execute()
{
     double aTime = cMMVII_Appli::CurrentAppli().SecFromT0();

     if (aTime < (mLastime+mPeriod))
        return false;

     mLastime += mPeriod;

     return true;
}


/***********************************/
/*                                 */
/*          cAutoTimerSegm         */
/*                                 */
/***********************************/

cAutoTimerSegm::cAutoTimerSegm(cTimerSegm & aTS,const tIndTS & anInd) :
   mTS        (&aTS),
   mSaveInd   (aTS.CurThreadState().mLastIndex),
   mBeginTime (aTS.mAppli->SecFromT0())
{
   mTS->SetIndex(anInd);
}

cAutoTimerSegm::cAutoTimerSegm(const tIndTS & anInd) :
    cAutoTimerSegm(GlobAppTS(),anInd)
{
}

cAutoTimerSegm::cAutoTimerSegm(cTimerSegm * aTS,const tIndTS & anInd) :
        mTS (aTS)
{
   if (aTS)
   {
      mSaveInd   =  aTS->CurThreadState().mLastIndex;
      mBeginTime =  aTS->mAppli->SecFromT0();
      mTS->SetIndex(anInd);
   }
}

cAutoTimerSegm::~cAutoTimerSegm()
{
   if (mTS)
   {
      // LIFO nesting guarantees this thread's open segment is the one we opened
      cTimerSegm::cThreadState & aTh = mTS->CurThreadState();
      tIndTS aClosedInd = aTh.mLastIndex;
      mTS->SetIndex(mSaveInd);
      aTh.mInclusive[aClosedInd] += mTS->mAppli->SecFromT0() - mBeginTime;
   }
}

/***********************************/
/*                                 */
/*          cTimerSegm             */
/*                                 */
/***********************************/

static const std::string DefTime("OTHERS");
static std::atomic<tU_INT8> TheNextTimerSegmId(1);

bool cTimerSegm::TheDefaultShowPerThread = false;

cTimerSegm::cThreadState::cThreadState(const tIndTS & aFirstIndex,double aT0) :
   mLastIndex     (aFirstIndex),
   mCurBeginTime  (aT0)
{
}

cTimerSegm::cTimerSegm(cMMVII_Ap_CPU * anAppli) :
   mId            (TheNextTimerSegmId.fetch_add(1)),
   mAppli         (anAppli),
   mShowAtDel     (true),
   mShowPerThread (TheDefaultShowPerThread)
{
}

// Per-thread cache, keyed by mId (not "this") so a stale entry can't survive a
// cTimerSegm destroyed and reallocated at the same address. State itself lives in
// mThreadStates, owned by the cTimerSegm, so it outlives the thread that filled it.
cTimerSegm::cThreadState & cTimerSegm::CurThreadState() const
{
   thread_local std::unordered_map<tU_INT8,cThreadState*> TheCache;

   auto anIt = TheCache.find(mId);
   if (anIt != TheCache.end())
      return *(anIt->second);

   std::lock_guard<std::mutex> aLock(mMutex);
   mThreadStates.emplace_back(DefTime,mAppli->SecFromT0());
   cThreadState * aRes = &mThreadStates.back();
   TheCache.emplace(mId,aRes);
   return *aRes;
}

double cTimerSegm::CurBeginTime() const {return  CurThreadState().mCurBeginTime;}

void  cTimerSegm::SetNoShowAtDel()
{
    mShowAtDel = false;
}

void  cTimerSegm::SetShowPerThread(bool aShow)
{
    mShowPerThread = aShow;
}

void  cTimerSegm::SetDefaultShowPerThread(bool aShow)
{
    TheDefaultShowPerThread = aShow;
}

cTimerSegm::~cTimerSegm()
{
   if (mShowAtDel && (Times().size() >=2)) // is something was added
     Show();
}

cTimerSegm & GlobAppTS()
{
   return cMMVII_Appli::CurrentAppli().TimeSegm();
}

tTableIndTS  cTimerSegm::Times() const
{
   std::lock_guard<std::mutex> aLock(mMutex);
   tTableIndTS aRes;
   for (const auto & aTh : mThreadStates)
       for (const auto & aPair : aTh.mExclusive)
           aRes[aPair.first] += aPair.second;
   return aRes;
}

tTableIndTS  cTimerSegm::TimesInclusive() const
{
   std::lock_guard<std::mutex> aLock(mMutex);
   tTableIndTS aRes;
   for (const auto & aTh : mThreadStates)
       for (const auto & aPair : aTh.mInclusive)
           aRes[aPair.first] += aPair.second;
   return aRes;
}

void cTimerSegm::SetIndex(const tIndTS & aInd)
{
   cThreadState & aTh = CurThreadState();
   double aCurTime =  mAppli->SecFromT0();
   aTh.mExclusive[aTh.mLastIndex] += aCurTime-aTh.mCurBeginTime;
   aTh.mCurBeginTime = aCurTime;
   aTh.mLastIndex = aInd;
}

// Names of aExcl, by decreasing duration -- systematic, no option
static std::vector<tIndTS> NamesByDecreasingDuration(const tTableIndTS & aExcl)
{
   std::vector<tIndTS> aRes;
   aRes.reserve(aExcl.size());
   for (const auto & aPair : aExcl)
       aRes.push_back(aPair.first);
   std::stable_sort
   (
       aRes.begin(),aRes.end(),
       [&aExcl](const tIndTS & aN1,const tIndTS & aN2) {return aExcl.at(aN1) > aExcl.at(aN2);}
   );
   return aRes;
}

void cTimerSegm::Show()
{
   {
      cAutoTimerSegm aATS(*this,DefTime);
   }

   tTableIndTS aExcl = Times();
   tTableIndTS anIncl = TimesInclusive();

   double aSom = 0.0;
   StdOut()  <<  Color::title << " ========== TIMING ===========" << Color::end << std::endl;
   for (const auto & aName : NamesByDecreasingDuration(aExcl))
   {
       double aDur = aExcl.at(aName);
       aSom += aDur;
       StdOut() << " * "  << FixDigToStr(aDur,4,4) << " : " << Color::command << aName << Color::end;
       auto anItIncl = anIncl.find(aName);
       // show inclusive time only when something ran nested inside this segment
       if ((anItIncl != anIncl.end()) && (anItIncl->second > aDur + 1e-3))
           StdOut() << Color::descr << "  (incl. nested: " << FixDigToStr(anItIncl->second,4,4) << ")" << Color::end;
       StdOut() << std::endl;
   }

   double aElapsed = mAppli->SecFromT0();
   StdOut() << " *** Total sum: " << aSom  <<  " Total ellapsed: " << aElapsed;
   // sum > elapsed is expected with concurrent threads, not an inconsistency
   if (aSom > aElapsed + 1e-3)
       StdOut() << Color::warning << "  (sum>elapsed : segments were measured concurrently by several threads, parallelism ratio="
                 << FixDigToStr(aSom/std::max(aElapsed,1e-9),3,2) << ")" << Color::end;
   StdOut() << std::endl;

   // per-thread breakdown, opt-in via SetShowPerThread
   std::lock_guard<std::mutex> aLock(mMutex);
   if (mShowPerThread && (mThreadStates.size() >= 2))
   {
       StdOut() << Color::title << " ---------- per thread ----------" << Color::end << std::endl;
       int aKTh = 0;
       for (const auto & aTh : mThreadStates)
       {
           StdOut() << Color::title << " -- Thread " << aKTh << " --" << Color::end << std::endl;
           for (const auto & aName : NamesByDecreasingDuration(aTh.mExclusive))
           {
               double aDur = aTh.mExclusive.at(aName);
               StdOut() << "    * " << FixDigToStr(aDur,4,4) << " : " << Color::command << aName << Color::end;
               auto anItIncl = aTh.mInclusive.find(aName);
               if ((anItIncl != aTh.mInclusive.end()) && (anItIncl->second > aDur + 1e-3))
                   StdOut() << Color::descr << "  (incl. nested: " << FixDigToStr(anItIncl->second,4,4) << ")" << Color::end;
               StdOut() << std::endl;
           }
           aKTh++;
       }
   }
}


};

