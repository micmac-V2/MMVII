#ifndef MMVII_CODETIMING_H
#define MMVII_CODETIMING_H

#include <chrono>
#include <string>
#include <iostream>
#include <iomanip>

namespace MMVII {
// ---------------------------------------------------------------------------
// cCodeTimer
//
// Measures the *cumulative* wall-clock time spent inside one or more
// instrumented code regions.  Each call to start()/stop() (or one
// ScopedTimer guard) adds a sample to the running total.  The timer can
// be paused and resumed any number of times — only the time between
// start() and stop() is counted.
//
// Clock choice — std::chrono::steady_clock:
//   • Monotonic: never goes backwards (safe across NTP adjustments, DST).
//   • Maps to RDTSC on x86 Linux/Windows/macOS → ~10–50 ns resolution.
//   • Do NOT use system_clock here; it can jump on clock corrections.
//
// Duration precision — nanoseconds:
//   Accumulation is done in nanoseconds to avoid rounding on every
//   stop().  Conversion to ms / µs only happens at reporting time.
// ---------------------------------------------------------------------------
class cCodeTimer {
public:
    using Clock    = std::chrono::steady_clock;
    using Duration = std::chrono::nanoseconds;

    explicit cCodeTimer(std::string name) : name_(std::move(name)) {}

    // Mark the beginning of a timed region.
    // Calling start() twice without an intervening stop() discards the
    // first timestamp — guard against this with ScopedTimer.
    void start() {
        t0_ = Clock::now();
    }

    // Mark the end of a timed region and add the elapsed time to the
    // running total.  The delta is computed as (now - t0_), so the
    // overhead of stop() itself is not included.
    void stop() {
        total_ += Clock::now() - t0_;
        ++calls_;
    }

    // Reset the accumulated time and call count to zero without
    // destroying the timer object (useful for per-frame stats).
    void reset() { total_ = Duration::zero(); calls_ = 0; }

    // --- Accessors -----------------------------------------------------------

    // Raw accumulated duration (nanosecond precision).
    Duration  total()    const { return total_; }

    // Number of start()/stop() pairs recorded so far.
    long long calls()    const { return calls_; }

    // Total time in milliseconds (double for sub-ms precision).
    double    total_ms() const {
        return std::chrono::duration<double, std::milli>(total_).count();
    }

    // Average time per call in microseconds.  Returns 0 if never called.
    double    avg_us()   const {
        if (calls_ == 0) return 0.0;
        return std::chrono::duration<double, std::micro>(total_).count() / calls_;
    }

    // Print a one-line summary to stdout.
    void report() const {
        std::cout << "[Timer] " << name_
                  << "  total="  << std::fixed << std::setprecision(3) << total_ms() << " ms"
                  << "  calls=" << calls_
                  << "  avg="   << std::setprecision(2) << avg_us()   << " µs\n";
    }

private:
    std::string       name_;
    Clock::time_point t0_{};    // timestamp of last start()
    Duration          total_{}; // accumulated wall time
    long long         calls_{}; // number of completed intervals
};

// ---------------------------------------------------------------------------
// cScopedCodeTimer  (RAII wrapper)
//
// Calls start() on construction and stop() on destruction, even if the
// enclosed scope exits via an exception.  This is the safest way to
// instrument a block: you cannot forget to call stop(), and early
// returns are handled automatically.
//
// Usage:
//   { cScopedCodeTimer g(myTimer);  do_work(); }  // stop() called here
// ---------------------------------------------------------------------------
struct cScopedCodeTimer {
    cCodeTimer& t;
    explicit cScopedCodeTimer(cCodeTimer& t) : t(t) { t.start(); }
    ~cScopedCodeTimer()                       { t.stop(); }
    // Non-copyable: copying a guard would double-stop the timer.
    cScopedCodeTimer(const cScopedCodeTimer&) = delete;
    cScopedCodeTimer& operator=(const cScopedCodeTimer&) = delete;
};


/*
  Usage Example:

cCodeTimer tA("parsing"), tB("rendering");

for (int i = 0; i < 1000; ++i) {

    // RAII style: stop() fires automatically at the closing brace,
    // even if do_parse() throws.
    { cScopedCodeTimer g(tA);  do_parse(i); }

    // Explicit style: useful when start/stop don't align with a scope.
    tB.start();
    do_render(i);
    tB.stop();
}

tA.report();  // [Timer] parsing    total=12.847 ms  calls=1000  avg=12.85 µs
tB.report();  // [Timer] rendering  total=38.012 ms  calls=1000  avg=38.01 µs
*/

} // namespace MMVII


#endif // MMVII_CODETIMING_H
