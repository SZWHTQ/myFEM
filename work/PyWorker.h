// This file does not belong to any project, just for the PyWorker interface.
#pragma once
#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif
EXPORT double PyWorker(const double L, const double B, const double ksi,
                       const double a_b, const double loadValue,
                       const double matrixElasticModulus,
                       const double matrixPoissonRatio,
                       const double inclusionElasticModulus,
                       const double inclusionPoissonRatio,
                       const double meshSize, const double refinementFactor,
                       const int meshAlgorithm, const bool isSerendipity,
                       const bool convertToSquare, const bool writeInp,
                       const bool writeMsh, const bool runFltk,
                       const bool isPlaneStress, const bool verbose);
#ifdef __cplusplus
}
#endif
