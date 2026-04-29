// Heavily inspired by https://github.com/PacktPublishing/Cpp-High-Performance/blob/master/Chapter03/scoped_timer.cpp


#include <chrono>

#define USE_TIMER 1

#if USE_TIMER
#define MEASURE_FUNCTION() ScopedTimer timer{__func__}
#else
#define MEASURE_FUNCTION()
#endif

class cScopedCodeTimer {

public:
    using ClockType = std::chrono::steady_clock;

  cScopedCodeTimer(const char* func);

  cScopedCodeTimer(const cScopedCodeTimer&) = delete;
  cScopedCodeTimer(cScopedCodeTimer&&) = delete;
  auto operator=(const cScopedCodeTimer&) -> cScopedCodeTimer& = delete;
  auto operator=(cScopedCodeTimer&&) -> cScopedCodeTimer& = delete;

  ~cScopedCodeTimer();

private:
  const char* function_ = {};
  const ClockType::time_point start_ = {};
};

