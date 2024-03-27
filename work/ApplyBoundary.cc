#include <iostream>
#include <memory>

#include "ApplyBoundary.h"
#include "Boundary.h"

std::list<Load> apply_load(Mesh* mesh, double L, double B) {
    std::list<Load> loadCondition;
    for (auto&& element : mesh->Elements) {
        if (element->elementName == "Serendipity") {
            for (int i = 0; i < 3; ++i) {
                // Upper bounds
                if (std::fabs(element->nodes[i]->x - L / 2) < 1e-6 &&
                    std::fabs(element->nodes[i + 1]->x - L / 2) < 1e-6) {
                    if (element->nodes[i]->x < element->nodes[i + 1]->x) {
                        Load load;
                        load.nodes[0] = element->nodes[i + 1];
                        load.nodes[1] = element->nodes[i];
                        load.nodes[2] = element->nodes[i + 4];
                        load.normalForce[0] = 1.0;
                        load.normalForce[1] = 1.0;
                        load.normalForce[2] = 1.0;
                        load.thickness = element->thickness;
                        loadCondition.push_back(load);
                    } else {
                        Load load;
                        load.nodes[0] = element->nodes[i];
                        load.nodes[1] = element->nodes[i + 1];
                        load.nodes[2] = element->nodes[i + 4];
                        load.normalForce[0] = 1.0;
                        load.normalForce[1] = 1.0;
                        load.normalForce[2] = 1.0;
                        load.thickness = element->thickness;
                        loadCondition.push_back(load);
                    }
                }
                // Right bounds
                if (std::fabs(element->nodes[i]->y - B / 2) < 1e-6 &&
                    std::fabs(element->nodes[i + 1]->y - B / 2) < 1e-6) {
                    if (element->nodes[i]->y < element->nodes[i + 1]->y) {
                        Load load;
                        load.nodes[0] = element->nodes[i + 1];
                        load.nodes[1] = element->nodes[i];
                        load.nodes[2] = element->nodes[i + 4];
                        load.normalForce[0] = 1.0;
                        load.normalForce[1] = 1.0;
                        load.normalForce[2] = 1.0;
                        load.thickness = element->thickness;
                        loadCondition.push_back(load);
                    } else {
                        Load load;
                        load.nodes[0] = element->nodes[i];
                        load.nodes[1] = element->nodes[i + 1];
                        load.nodes[2] = element->nodes[i + 4];
                        load.normalForce[0] = 1.0;
                        load.normalForce[1] = 1.0;
                        load.normalForce[2] = 1.0;
                        load.thickness = element->thickness;
                        loadCondition.push_back(load);
                    }
                }
                // // Lower bounds
                // if (std::fabs(element->nodes[i]->x + L / 2) < 1e-6 &&
                //     std::fabs(element->nodes[i + 1]->x + L / 2) < 1e-6) {
                //     if (element->nodes[i]->x < element->nodes[i + 1]->x) {
                //         Load load;
                //         load.nodes[0] = element->nodes[i + 1];
                //         load.nodes[1] = element->nodes[i];
                //         load.nodes[2] = element->nodes[i + 4];
                //         load.normalForce[0] = -1.0;
                //         load.normalForce[1] = -1.0;
                //         load.normalForce[2] = -1.0;
                //         load.thickness = element->thickness;
                //         loadCondition.push_back(load);
                //     } else {
                //         Load load;
                //         load.nodes[0] = element->nodes[i];
                //         load.nodes[1] = element->nodes[i + 1];
                //         load.nodes[2] = element->nodes[i + 4];
                //         load.normalForce[0] = -1.0;
                //         load.normalForce[1] = -1.0;
                //         load.normalForce[2] = -1.0;
                //         load.thickness = element->thickness;
                //         loadCondition.push_back(load);
                //     }
                // }
                // // Left bounds
                // if (std::fabs(element->nodes[i]->y + B / 2) < 1e-6 &&
                //     std::fabs(element->nodes[i + 1]->y + B / 2) < 1e-6) {
                //     if (element->nodes[i]->y < element->nodes[i + 1]->y) {
                //         Load load;
                //         load.nodes[0] = element->nodes[i];
                //         load.nodes[1] = element->nodes[i + 1];
                //         load.nodes[2] = element->nodes[i + 4];
                //         load.normalForce[0] = -1.0;
                //         load.normalForce[1] = -1.0;
                //         load.normalForce[2] = -1.0;
                //         load.thickness = element->thickness;
                //         loadCondition.push_back(load);
                //     } else {
                //         Load load;
                //         load.nodes[0] = element->nodes[i + 1];
                //         load.nodes[1] = element->nodes[i];
                //         load.nodes[2] = element->nodes[i + 4];
                //         load.normalForce[0] = -1.0;
                //         load.normalForce[1] = -1.0;
                //         load.normalForce[2] = -1.0;
                //         load.thickness = element->thickness;
                //         loadCondition.push_back(load);
                //     }
                // }
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
        if (std::fabs(node.x - X) < 1e-6) {
            Boundary boundary;
            boundary.nodes = std::make_shared<Node>(node);
            boundary.fixed[0] = true;
            boundaryCondition.push_back(boundary);
        }
        if (std::fabs(node.y - Y) < 1e-6) {
            Boundary boundary;
            boundary.nodes = std::make_shared<Node>(node);
            boundary.fixed[1] = true;
            boundaryCondition.push_back(boundary);
        }
    }
    return boundaryCondition;
}