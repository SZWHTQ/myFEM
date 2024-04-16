#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT __attribute__((visibility("default")))
#endif

// extern "C" {
// EXPORT double work(char* tomlFilePath, double inclusionElasticModulus,
//                    double ksi, double a_b, bool verbose);
// }
