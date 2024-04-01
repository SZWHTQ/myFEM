#pragma once
// #include <time.h>
#include <chrono>

// Millisecond timer
class Timer {
   public:
    Timer() : start(std::chrono::high_resolution_clock::now()) {}
    void reset() { start = std::chrono::high_resolution_clock::now(); }
    double elapsed() {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(end - start).count();
    }

   private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
};