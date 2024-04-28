#pragma once

#include <chrono>
#include <ostream>

// Millisecond timer
class Timer {
   public:
    Timer() : start(std::chrono::high_resolution_clock::now()) {}
    ~Timer() = default;
    void reset() { start = std::chrono::high_resolution_clock::now(); }
    [[nodiscard]] double elapsed() const {
        auto now = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(now - start).count();
    }

    friend std::ostream& operator<<(std::ostream& os, const Timer& timer);

   private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
};