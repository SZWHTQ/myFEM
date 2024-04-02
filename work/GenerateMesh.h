#pragma once

#include <vector>

#define WRITE_INP 1
#define WRITE_MSH 0

int generate_mesh(std::vector<double>& nodeCoord,
                  std::vector<size_t>& elementNodeTags, double L, double B,
                  double a, double b, double lc, double refinementFactor = 1,
                  bool isSerendipity = true, int meshAlgorithm = 8);