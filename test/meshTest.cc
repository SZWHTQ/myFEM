#include <cmath>
#include <iostream>

#include "GenerateMesh.h"

int main() {
    std::vector<double> nodeCoord;
    std::vector<size_t> elementNodeTags;
    double L = 2;
    double B = 0.9;
    double a = 0.3;
    double b = 0.1;
    double lc = 0.2;
    int val = generate_mesh(nodeCoord, elementNodeTags, L, B, a, b, lc);
    if (val != 0) {
        return -1;
    }
    for (size_t i = 0; i < elementNodeTags.size(); ++i) {
        size_t tag = elementNodeTags[i];

        if (std::fabs(nodeCoord[3 * (tag - 1)] - L / 2) < 1e-6 * lc) {
            std::cout << "Element " << i / 8 + 1
                      << " is on the right boundary ";
            std::cout << "Node " << i % 8 + 1 << ", " << tag << std::endl;
        }
        if (std::fabs(nodeCoord[3 * (tag - 1)] - (-L / 2)) < 1e-6 * lc) {
            std::cout << "Element " << i / 8 + 1 << " is on the left boundary ";
            std::cout << "Node " << i % 8 + 1 << ", " << tag << std::endl;
        }
        if (std::fabs(nodeCoord[3 * (tag - 1) + 1] - B / 2) < 1e-6 * lc) {
            std::cout << "Element " << i / 8 + 1 << " is on the top boundary ";
            std::cout << "Node " << i % 8 + 1 << ", " << tag << std::endl;
        }
        if (std::fabs(nodeCoord[3 * (tag - 1) + 1] - (-B / 2)) < 1e-6 * lc) {
            std::cout << "Element " << i / 8 + 1
                      << " is on the bottom boundary ";
            std::cout << "Node " << i % 8 + 1 << ", " << tag << std::endl;
        }
    }
}