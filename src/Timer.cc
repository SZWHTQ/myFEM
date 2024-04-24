#include <cmath>
#include <iomanip>

#include "Timer.h"

std::ostream& operator<<(std::ostream& os, const Timer& timer) {
    os << std::setprecision(3);
    double milliSeconds = timer.elapsed();
    if (milliSeconds < 1e3) {
        os << milliSeconds << "ms";
    } else {
        double seconds = milliSeconds / 1e3;
        if (seconds < 60.0) {
            os << seconds << "s";
        } else {
            unsigned int minutes = (unsigned int)seconds / 60;
            seconds = fmod(seconds, 60.0);
            if (minutes < 60) {
                os << minutes << "m " << seconds << "s";
            } else {
                unsigned int hours = minutes / 60;
                minutes %= 60;
                os << hours << "h " << minutes << "m " << seconds << "s";
            }
        }
    }
    os << std::setprecision(-1);
    return os;
}