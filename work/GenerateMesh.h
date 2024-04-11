#pragma once

#include <vector>

#define WRITE_INP 1
#define WRITE_MSH 1

#define VERTEX_TAG 1
#define BOUNDARY_TAG 2
#define MATRIX_SURFACE_TAG 3
#define INCLUSION_SURFACE_TAG 4
#define INTERFACE_TAG 5

int generate_mesh(std::vector<double>& nodeCoord,
                  std::vector<size_t>& elementNodeTags,
                  std::vector<size_t>& boundaryNodeTags, double L, double B,
                  double a, double b, double lc, double refinementFactor = 1,
                  bool isSerendipity = true, int meshAlgorithm = 8);