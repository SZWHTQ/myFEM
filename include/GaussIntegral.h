#include <cmath>
#include <vector>

namespace GaussIntegral {
template <int N>
struct GaussData;

template <>
struct GaussData<2> {
    static inline const std::vector<double> weights = {1.0, 1.0};
    static inline const std::vector<double> abscissas = {-1.0 / std::sqrt(3.0),
                                                         1.0 / std::sqrt(3.0)};
};

template <>
struct GaussData<3> {
    static inline const std::vector<double> weights = {5.0 / 9.0, 8.0 / 9.0,
                                                       5.0 / 9.0};
    static inline const std::vector<double> abscissas = {
        -std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0)};
};

template <>
struct GaussData<4> {
    static inline const std::vector<double> weights = {
        (18.0 - std::sqrt(30.0)) / 36.0, (18.0 + std::sqrt(30.0)) / 36.0,
        (18.0 + std::sqrt(30.0)) / 36.0, (18.0 - std::sqrt(30.0)) / 36.0};
    static inline const std::vector<double> abscissas = {
        -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
        -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
        std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
        std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0))};
};

template <>
struct GaussData<5> {
    static inline const std::vector<double> weights = {
        (322.0 - 13.0 * std::sqrt(70.0)) / 900.0,
        (322.0 + 13.0 * std::sqrt(70.0)) / 900.0, 128.0 / 225.0,
        (322.0 + 13.0 * std::sqrt(70.0)) / 900.0,
        (322.0 - 13.0 * std::sqrt(70.0)) / 900.0};
    static inline const std::vector<double> abscissas = {
        -std::sqrt(5.0 + 2.0 * std::sqrt(10.0 / 7.0)) / 3.0,
        -std::sqrt(5.0 - 2.0 * std::sqrt(10.0 / 7.0)) / 3.0, 0.0,
        std::sqrt(5.0 - 2.0 * std::sqrt(10.0 / 7.0)) / 3.0,
        std::sqrt(5.0 + 2.0 * std::sqrt(10.0 / 7.0)) / 3.0};
};

template <>
struct GaussData<6> {
    static inline const std::vector<double> weights = {
        0.171324492379170345040296142172732893026705587812571628056,
        0.360761573048138607569833513837716111661694098041216356020,
        0.467913934572691047389870343989550994811655505647615320970,
        0.467913934572691047389870343989550994811655505647615320970,
        0.360761573048138607569833513837716111661694098041216356020,
        0.171324492379170345040296142172732893026705587812571628056};
    static inline const std::vector<double> abscissas = {
        -0.932469514203152027812301554493994609134765545647383745404,
        -0.661209386466264513661399595019905347006448594769820976748,
        -0.238619186083196908630501721680711935418610043146485340484,
        0.238619186083196908630501721680711935418610043146485340484,
        0.661209386466264513661399595019905347006448594769820976748,
        0.932469514203152027812301554493994609134765545647383745404};
};

}  // namespace GaussIntegral