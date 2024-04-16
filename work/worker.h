#pragma once
#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif
EXPORT double worker(char* tomlFilePath, double inclusionElasticModulus,
                   double ksi, double a_b, bool verbose);
#ifdef __cplusplus
}
#endif
