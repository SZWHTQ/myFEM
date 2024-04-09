#include <gmsh.h>

#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>

#include "GetStrainEnergyChange.h"
#include "Material.h"

double getDistance(std::shared_ptr<Node> node1, std::shared_ptr<Node> node2) {
    return std::sqrt(std::pow(node1->x - node2->x, 2) +
                     std::pow(node1->y - node2->y, 2) +
                     std::pow(node1->z - node2->z, 2));
}

Eigen::Vector2d getNormal(std::shared_ptr<Node> node1,
                          std::shared_ptr<Node> node2) {
    Eigen::Vector2d normal(node2->y - node1->y, node1->x - node2->x);
    normal.normalize();
    return normal;
}

Eigen::Vector2d getTraction(std::shared_ptr<Node> node1,
                            std::shared_ptr<Node> node2,
                            std::shared_ptr<Node> node3) {
    Eigen::Vector2d traction;
    Eigen::Matrix2d stress;
    stress(0, 0) = (node1->Stress[0] + node2->Stress[0] + node3->Stress[0]) / 3;
    stress(0, 1) = (node1->Stress[2] + node2->Stress[2] + node3->Stress[2]) / 3;
    stress(1, 0) = (node1->Stress[2] + node2->Stress[2] + node3->Stress[2]) / 3;
    stress(1, 1) = (node1->Stress[1] + node2->Stress[1] + node3->Stress[1]) / 3;
    auto&& normal = getNormal(node1, node2);
    traction(0) = stress(0, 0) * normal(0) + stress(0, 1) * normal(1);
    traction(1) = stress(1, 0) * normal(0) + stress(1, 1) * normal(1);

    return traction;
}

double getStressSum(std::shared_ptr<Node> node1, std::shared_ptr<Node> node2,
                    std::shared_ptr<Node> node3) {
    double stressSum = 0;
    stressSum += node1->Stress[0] + node1->Stress[1];
    stressSum += node2->Stress[0] + node2->Stress[1];
    stressSum += node3->Stress[0] + node3->Stress[1];
    stressSum /= 3;
    return stressSum;
}

Eigen::Vector2d getDisplacement(std::shared_ptr<Node> node1,
                                std::shared_ptr<Node> node2,
                                std::shared_ptr<Node> node3) {
    Eigen::Vector2d displacement;
    displacement(0) = (node1->Displacement[0] + node2->Displacement[0] +
                       node3->Displacement[0]) /
                      3;
    displacement(1) = (node1->Displacement[1] + node2->Displacement[1] +
                       node3->Displacement[1]) /
                      3;
    return displacement;
}

double getStrainEnergyChange(Mesh* mesh, Mesh* meshNoInclusion,
                             Material* matrix, Material* inclusion) {
    double deltaU = 0;
    double lambda1 =
        1 -
        ((1 + matrix->getPoissonRatio()) * inclusion->getElasticModulus()) /
            ((1 + inclusion->getPoissonRatio()) / matrix->getElasticModulus());
    double lambda2 =
        ((inclusion->getPoissonRatio() - matrix->getPoissonRatio()) *
         (1 + matrix->getPoissonRatio()) * inclusion->getElasticModulus()) /
        ((1 + inclusion->getPoissonRatio()) *
         (1 - 2 * inclusion->getPoissonRatio()) * matrix->getElasticModulus());
    // double lambda1 = 1;
    // double lambda2 = 1;

    for (auto& element : meshNoInclusion->Elements) {
        size_t elementIndex = element->getIndex();
        if (mesh->Elements[elementIndex]->material->getIndex() == 1) {
            for (int i = 0; i < 4; ++i) {
                int j = (i + 1) % 4;
                auto& n1 = element->nodes[i];
                auto& n2 = element->nodes[j];
                auto& n3 = element->nodes[i + 4];
                if (n1->isBoundary && n2->isBoundary) {
                    auto&& distance = getDistance(n1, n2);
                    auto&& traction = getTraction(n1, n2, n3);
                    auto&& stressSum = getStressSum(n1, n2, n3);
                    auto&& normal = getNormal(n1, n2);
                    auto&& displacement =
                        getDisplacement(mesh->Nodes[n1->getIndex()],
                                        mesh->Nodes[n2->getIndex()],
                                        mesh->Nodes[n3->getIndex()]);
                    deltaU +=
                        0.5 *
                        (lambda1 * traction - lambda2 * stressSum * normal)
                            .dot(displacement) *
                        distance;
                }
            }
        }
    }
    return deltaU;
}