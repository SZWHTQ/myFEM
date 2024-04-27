#include <iostream>
#include <memory>

#include "ApplyBoundary.h"
#include "Boundary.h"
#include "Element.h"

std::list<Load> apply_load(Mesh* mesh, double L, double B, double value) {
    std::list<Load> loadCondition;
    for (auto&& element : mesh->Elements) {
        if (element->getElementName() == "Serendipity") {
            for (int i = 0; i < 4; ++i) {
                int j = (i + 1) % 4;
                // Right bound
                if (std::fabs(element->nodes[i]->x - L / 2) < 1e-6 &&
                    std::fabs(element->nodes[j]->x - L / 2) < 1e-6) {
                    if (element->nodes[i]->y > element->nodes[j]->y) {
                        Load load;
                        load.nodes[0] = element->nodes[j];
                        load.nodes[1] = element->nodes[i];
                        load.nodes[2] = element->nodes[i + 4];
                        load.normalForce[0] = value;
                        load.normalForce[1] = value;
                        load.normalForce[2] = value;
                        load.thickness = element->thickness;
                        loadCondition.push_back(load);
                    } else {
                        Load load;
                        load.nodes[0] = element->nodes[i];
                        load.nodes[1] = element->nodes[j];
                        load.nodes[2] = element->nodes[i + 4];
                        load.normalForce[0] = value;
                        load.normalForce[1] = value;
                        load.normalForce[2] = value;
                        load.thickness = element->thickness;
                        loadCondition.push_back(load);
                    }
                }
                // Upper bound
                if (std::fabs(element->nodes[i]->y - B / 2) < 1e-6 &&
                    std::fabs(element->nodes[j]->y - B / 2) < 1e-6) {
                    if (element->nodes[i]->x < element->nodes[j]->x) {
                        Load load;
                        load.nodes[0] = element->nodes[j];
                        load.nodes[1] = element->nodes[i];
                        load.nodes[2] = element->nodes[i + 4];
                        load.normalForce[0] = value;
                        load.normalForce[1] = value;
                        load.normalForce[2] = value;
                        load.thickness = element->thickness;
                        loadCondition.push_back(load);
                    } else {
                        Load load;
                        load.nodes[0] = element->nodes[i];
                        load.nodes[1] = element->nodes[j];
                        load.nodes[2] = element->nodes[i + 4];
                        load.normalForce[0] = value;
                        load.normalForce[1] = value;
                        load.normalForce[2] = value;
                        load.thickness = element->thickness;
                        loadCondition.push_back(load);
                    }
                }
            }
        } else {
            std::cout << "Element type not supported\n";
        }
    }
    return loadCondition;
}

std::list<Boundary> apply_boundary(Mesh* mesh, double X, double Y) {
    std::list<Boundary> boundaryCondition;
    for (auto&& node : mesh->Nodes) {
        if (std::fabs(node->x - X) < 1e-6) {
            Boundary boundary;
            boundary.node = node;
            boundary.fixed[0] = true;
            boundaryCondition.push_back(boundary);
        }
        if (std::fabs(node->y - Y) < 1e-6) {
            Boundary boundary;
            boundary.node = node;
            boundary.fixed[1] = true;
            boundaryCondition.push_back(boundary);
        }
    }
    return boundaryCondition;
}