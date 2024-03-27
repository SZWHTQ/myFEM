#pragma once
#include <stdexcept>
#include <vector>

class GaussIntegral {
   public:
    struct GaussData {
        std::vector<double> weights;
        std::vector<double> abscissas;
    };

   private:
    // Private constructor to prevent instantiation
    GaussIntegral() {}

    // GaussData for different point numbers
    static const GaussData twoPoint;
    static const GaussData threePoint;
    static const GaussData fourPoint;
    static const GaussData fivePoint;
    static const GaussData sixPoint;

   public:
    // Public static method to access GaussData

    static const GaussData& getGaussData(int gaussPointNum) {
        switch (gaussPointNum) {
            case 2:
                return twoPoint;
            case 3:
                return threePoint;
            case 4:
                return fourPoint;
            case 5:
                return fivePoint;
            case 6:
                return sixPoint;
            default:
                throw std::invalid_argument("Unsupported number of points");
        }
    }
};
