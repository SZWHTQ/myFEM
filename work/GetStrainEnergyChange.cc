#include <gmsh.h>

#include <memory>
#include <vector>

#include "GetStrainEnergyChange.h"
#include "Material.h"

static double getDistance(std::shared_ptr<Node> node1,
                          std::shared_ptr<Node> node2) {
    return std::sqrt(std::pow(node1->x - node2->x, 2) +
                     std::pow(node1->y - node2->y, 2) +
                     std::pow(node1->z - node2->z, 2));
}

static Eigen::Vector2d getNormal(std::shared_ptr<Node> node1,
                                 std::shared_ptr<Node> node2) {
    Eigen::Vector2d normal(node2->y - node1->y, node1->x - node2->x);
    normal.normalize();
    return normal;
}

static Eigen::Vector2d getTraction(std::shared_ptr<Node> node1,
                                   std::shared_ptr<Node> node2) {
    Eigen::Vector2d traction;
    auto&& normal = getNormal(node1, node2);
    traction(0) = 1 * normal(0);
    traction(1) = 1 * normal(1);

    return traction;
}

inline static double getStressSum() { return 2; }

static Eigen::Vector2d getDisplacement(std::shared_ptr<Node> node1,
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

double getStrainEnergyChange(Mesh* mesh, Material* matrix, Material* inclusion,
                             bool isPlaneStress) {
    double deltaU = 0;
    double lambda1 =
        1 -
        ((1 + matrix->getPoissonRatio()) * inclusion->getElasticModulus()) /
            ((1 + inclusion->getPoissonRatio()) / matrix->getElasticModulus());

    double temp = 1;
    if (!isPlaneStress) {
        temp = 1 + matrix->getPoissonRatio();
    }

    double lambda2 =
        ((inclusion->getPoissonRatio() - matrix->getPoissonRatio()) * temp *
         inclusion->getElasticModulus()) /
        ((1 + inclusion->getPoissonRatio()) *
         (1 - 2 * inclusion->getPoissonRatio()) * matrix->getElasticModulus());
    // double lambda1 = 1;
    // double lambda2 = 1;

    for (auto& element : mesh->Elements) {
        size_t elementIndex = element->getIndex();
        if (element->material->getIndex() == 1) {
            for (int i = 0; i < 4; ++i) {
                int j = (i + 1) % 4;
                auto& n1 = element->nodes[i];
                auto& n2 = element->nodes[j];
                auto& n3 = element->nodes[i + 4];
                if (n1->isBoundary && n2->isBoundary) {
                    auto&& distance = getDistance(n1, n2);
                    auto&& traction = getTraction(n1, n2);
                    auto&& stressSum = getStressSum();
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