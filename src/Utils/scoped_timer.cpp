#include <iostream>
#include "scoped_timer.h"

cScopedCodeTimer::cScopedCodeTimer(const char* func)
    : function_{func}, start_{ClockType::now()}
{
}

cScopedCodeTimer::~cScopedCodeTimer()
{
    using namespace std::chrono;
    auto stop = ClockType::now();
    auto duration = (stop - start_);
    auto ms = duration_cast<milliseconds>(duration).count();
    std::cout << ms << " ms " << function_ <<  '\n';
}

