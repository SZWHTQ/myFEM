#include <Eigen/SparseLU>
// #include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

#include "Element.h"
#include "Mesh.h"
#include "Serendipity.h"

Mesh::Mesh(MeshType type, std::vector<double> nodeCoord,
           std::vector<size_t> elementNodeTags, bool planeStress_)
    : meshType(type), planeStress(planeStress_) {
    for (size_t i = 0; i < nodeCoord.size(); i += 3) {
        Nodes.push_back(
            Node(i / 3, nodeCoord[i], nodeCoord[i + 1], nodeCoord[i + 2]));
    }
    // size_t minTag =
    //     *std::min_element(elementNodeTags.begin(), elementNodeTags.end());
    // std::cout << "Min tag: " << minTag << std::endl;
    switch (meshType) {
        case MeshType::serendipity:
            for (size_t i = 0; i < elementNodeTags.size();
                 i += Serendipity::nodeNum) {
                std::vector<std::shared_ptr<Node>> singleElementNodes;
                for (size_t j = 0; j < Serendipity::nodeNum; ++j) {
                    singleElementNodes.push_back(std::make_shared<Node>(
                        Nodes[elementNodeTags[i + j] - 1]));
                }
                Element* element = new Serendipity(
                    i / Serendipity::nodeNum, singleElementNodes, planeStress);
                Elements.push_back(element);
            }
            break;
        default:
            std::cerr << "Mesh type not implemented" << std::endl;
            break;
    }
}

const std::vector<double> Mesh::equivalentForce(Load* load) {
    std::vector<double> equivalentForce(6);
    Eigen::VectorXd X(6), Y(6);
    auto&& n = load->nodes;
    X(0) = -10 * n[0]->x - 2 * n[1]->x + 12 * n[2]->x;
    X(1) = n[0]->x - n[1]->x;
    X(2) = -6 * n[0]->x - 2 * n[1]->x + 8 * n[2]->x;
    X(3) = 2 * n[0]->x + 10 * n[1]->x - 12 * n[2]->x;
    X(4) = 2 * n[0]->x + 6 * n[1]->x - 8 * n[2]->x;
    X(5) = -16 * n[0]->x + 16 * n[1]->x;
    Y(0) = -10 * n[0]->y - 2 * n[1]->y + 12 * n[2]->y;
    Y(1) = n[0]->y - n[1]->y;
    Y(2) = -6 * n[0]->y - 2 * n[1]->y + 8 * n[2]->y;
    Y(3) = 2 * n[0]->y + 10 * n[1]->y - 12 * n[2]->y;
    Y(4) = 2 * n[0]->y + 6 * n[1]->y - 8 * n[2]->y;
    Y(5) = -16 * n[0]->y + 16 * n[1]->y;

    equivalentForce[0] =
        (X(0) * load->shearForce[0] + X(1) * load->shearForce[1] +
         X(2) * load->shearForce[2] + Y(0) * load->normalForce[0] +
         Y(1) * load->normalForce[1] + Y(2) * load->normalForce[2]) *
        load->thickness / 30;
    equivalentForce[1] =
        (Y(0) * load->shearForce[0] + Y(1) * load->shearForce[1] +
         Y(2) * load->shearForce[2] - X(0) * load->normalForce[0] -
         X(1) * load->normalForce[1] - X(2) * load->normalForce[2]) *
        load->thickness / 30;
    equivalentForce[2] =
        (X(1) * load->shearForce[0] + X(3) * load->shearForce[1] +
         X(4) * load->shearForce[2] + Y(1) * load->normalForce[0] +
         Y(3) * load->normalForce[1] + Y(4) * load->normalForce[2]) *
        load->thickness / 30;
    equivalentForce[3] =
        (Y(1) * load->shearForce[0] + Y(3) * load->shearForce[1] +
         Y(4) * load->shearForce[2] - X(1) * load->normalForce[0] -
         X(3) * load->normalForce[1] - X(4) * load->normalForce[2]) *
        load->thickness / 30;
    equivalentForce[4] =
        (X(2) * load->shearForce[0] + X(4) * load->shearForce[1] +
         X(5) * load->shearForce[2] + Y(2) * load->normalForce[0] +
         Y(4) * load->normalForce[1] + Y(5) * load->normalForce[2]) *
        load->thickness / 30;
    equivalentForce[5] =
        (Y(2) * load->shearForce[0] + Y(4) * load->shearForce[1] +
         Y(5) * load->shearForce[2] - X(2) * load->normalForce[0] -
         X(4) * load->normalForce[1] - X(5) * load->normalForce[2]) *
        load->thickness / 30;
    return equivalentForce;
}

Eigen::MatrixXd Mesh::assembleStiffnessMatrix() {
    Eigen::MatrixXd globalStiffnessMatrix(Nodes.size() * 2, Nodes.size() * 2);
    globalStiffnessMatrix.setZero();
    for (auto&& element : Elements) {
        auto&& ke = element->stiffnessMatrix();
        for (size_t i = 0; i < Serendipity::nodeNum; ++i) {
            for (size_t j = 0; j < Serendipity::nodeNum; ++j) {
                globalStiffnessMatrix(2 * element->nodes[i]->index,
                                      2 * element->nodes[j]->index) +=
                    ke(2 * i, 2 * j);
                globalStiffnessMatrix(2 * element->nodes[i]->index,
                                      2 * element->nodes[j]->index + 1) +=
                    ke(2 * i, 2 * j + 1);
                globalStiffnessMatrix(2 * element->nodes[i]->index + 1,
                                      2 * element->nodes[j]->index) +=
                    ke(2 * i + 1, 2 * j);
                globalStiffnessMatrix(2 * element->nodes[i]->index + 1,
                                      2 * element->nodes[j]->index + 1) +=
                    ke(2 * i + 1, 2 * j + 1);
            }
        }
    }
    return globalStiffnessMatrix;
}

int Mesh::Solve(std::list<Load>& loads, std::list<Boundary>& boundaries) {
    Force.resize(Nodes.size() * 2);
    Force.setZero();
    for (auto&& load : loads) {
        auto&& equivalentForce = Mesh::equivalentForce(&load);
        for (size_t i = 0; i < 3; ++i) {
            Force(2 * load.nodes[i]->index) += equivalentForce[2 * i];
            Force(2 * load.nodes[i]->index + 1) += equivalentForce[2 * i + 1];
        }
    }
    auto&& globalStiffnessMatrix = assembleStiffnessMatrix();
    // std::cout << "Global Stiffness Matrix:\n";
    // std::cout << globalStiffnessMatrix << std::endl;

    // Apply boundary conditions
    for (auto&& boundary : boundaries) {
        for (size_t i = 0; i < 2; ++i) {
            if (boundary.fixed[i]) {
                size_t j = 2 * boundary.node->index + i;
                globalStiffnessMatrix(j, j) *= 1e15;
                Force(j) = 0;
            }
        }
    }

    auto&& K = globalStiffnessMatrix.sparseView();
    auto&& F = Force.sparseView();
    Eigen::SparseVector<double> U;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(K);
    solver.factorize(K);
    if (solver.info() != Eigen::Success) {
        std::cerr << "LU decomposition failed" << std::endl;
        std::cerr << "Error code: " << solver.info() << std::endl;
        std::cerr << "Error message: " << solver.lastErrorMessage() << "\n";
        return -1;
    }
    U = solver.solve(F);

    for (size_t i = 0; i < Nodes.size(); ++i) {
        Nodes[i].displacement[0] = U.coeffRef(2 * i);
        Nodes[i].displacement[1] = U.coeffRef(2 * i + 1);
    }
    return 0;
}