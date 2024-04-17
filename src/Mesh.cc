#include <Eigen/SparseLU>
#include <iostream>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

#include "Element.h"
#include "Mesh.h"
#include "Serendipity.h"
#include "ThreadPool.h"
#include "Timer.h"

size_t Serendipity::nodeNum = 8;

Mesh::Mesh(MeshType type, std::vector<double> nodeCoord,
           std::vector<size_t> elementNodeTags,
           std::vector<size_t> boundaryNodeTags, bool planeStress_)
    : meshType(type), planeStress(planeStress_) {
    Nodes.reserve(nodeCoord.size() / 3);
    for (size_t i = 0; i < nodeCoord.size(); i += 3) {
        size_t id = i / 3;
        Nodes.push_back(std::make_shared<Node>(
            id, nodeCoord[i], nodeCoord[i + 1], nodeCoord[i + 2]));
    }

    for (int i = 0; i < boundaryNodeTags.size(); ++i) {
        Nodes[boundaryNodeTags[i]]->isBoundary = true;
    }

    switch (meshType) {
        case MeshType::serendipity:
            for (size_t i = 0; i < elementNodeTags.size();
                 i += Serendipity::nodeNum) {
                std::vector<std::shared_ptr<Node>> singleElementNodes;
                singleElementNodes.reserve(Serendipity::nodeNum);
                for (size_t j = 0; j < Serendipity::nodeNum; ++j) {
                    singleElementNodes.emplace_back(
                        Nodes[elementNodeTags[i + j] - 1]);
                }
                Element* element = new Serendipity(
                    i / Serendipity::nodeNum, singleElementNodes, planeStress);
                Elements.emplace_back(element);
            }
            break;
        default:
            std::cerr << "Mesh type not implemented" << std::endl;
            break;
    }
}

Mesh::Mesh(std::vector<double> nodeCoord,
           std::vector<std::pair<MeshType, std::vector<size_t>>>
               elementTypeAndNodeTags,
           std::vector<size_t> boundaryNodeTags, bool planeStress_)
    : planeStress(planeStress_) {
    Nodes.reserve(nodeCoord.size() / 3);
    for (size_t i = 0; i < nodeCoord.size(); i += 3) {
        Nodes.push_back(std::make_shared<Node>(
            i / 3, nodeCoord[i], nodeCoord[i + 1], nodeCoord[i + 2]));
    }

    for (int i = 0; i < boundaryNodeTags.size(); ++i) {
        Nodes[boundaryNodeTags[i] - 1]->isBoundary = true;
    }

    for (auto&& [type, elementNodeTags] : elementTypeAndNodeTags) {
        switch (type) {
            case MeshType::serendipity:
                for (size_t i = 0; i < elementNodeTags.size();
                     i += Serendipity::nodeNum) {
                    std::vector<std::shared_ptr<Node>> singleElementNodes;
                    singleElementNodes.reserve(Serendipity::nodeNum);
                    for (size_t j = 0; j < Serendipity::nodeNum; ++j) {
                        singleElementNodes.emplace_back(
                            Nodes[elementNodeTags[i + j] - 1]);
                    }
                    Element* element =
                        new Serendipity(i / Serendipity::nodeNum,
                                        singleElementNodes, planeStress);
                    Elements.emplace_back(element);
                }
                break;
            default:
                std::cerr << "Mesh type not implemented" << std::endl;
                break;
        }
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

    ThreadPool pool(std::thread::hardware_concurrency());
    for (auto element : Elements) {
        auto&& ke = element->getStiffnessMatrix();
        for (size_t i = 0; i < Serendipity::nodeNum; ++i) {
            for (size_t j = 0; j < Serendipity::nodeNum; ++j) {
                globalStiffnessMatrix(2 * element->nodes[i]->getIndex(),
                                      2 * element->nodes[j]->getIndex()) +=
                    ke(2 * i, 2 * j);
                globalStiffnessMatrix(2 * element->nodes[i]->getIndex(),
                                      2 * element->nodes[j]->getIndex() + 1) +=
                    ke(2 * i, 2 * j + 1);
                globalStiffnessMatrix(2 * element->nodes[i]->getIndex() + 1,
                                      2 * element->nodes[j]->getIndex()) +=
                    ke(2 * i + 1, 2 * j);
                globalStiffnessMatrix(2 * element->nodes[i]->getIndex() + 1,
                                      2 * element->nodes[j]->getIndex() + 1) +=
                    ke(2 * i + 1, 2 * j + 1);
            }
        }
    }

    return globalStiffnessMatrix;
}

Eigen::SparseMatrix<double> Mesh::sparseAssembleStiffnessMatrix() {
    Eigen::SparseMatrix<double> globalStiffnessMatrix(Nodes.size() * 2,
                                                      Nodes.size() * 2);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(Elements.size() * 4 * Serendipity::nodeNum *
                     Serendipity::nodeNum);
    for (auto element : Elements) {
        auto&& ke = element->getStiffnessMatrix();
        for (size_t i = 0; i < Serendipity::nodeNum; ++i) {
            for (size_t j = 0; j < Serendipity::nodeNum; ++j) {
                triplets.emplace_back(2 * element->nodes[i]->getIndex(),
                                      2 * element->nodes[j]->getIndex(),
                                      ke(2 * i, 2 * j));
                triplets.emplace_back(2 * element->nodes[i]->getIndex(),
                                      2 * element->nodes[j]->getIndex() + 1,
                                      ke(2 * i, 2 * j + 1));
                triplets.emplace_back(2 * element->nodes[i]->getIndex() + 1,
                                      2 * element->nodes[j]->getIndex(),
                                      ke(2 * i + 1, 2 * j));
                triplets.emplace_back(2 * element->nodes[i]->getIndex() + 1,
                                      2 * element->nodes[j]->getIndex() + 1,
                                      ke(2 * i + 1, 2 * j + 1));
            }
        }
    }

    globalStiffnessMatrix.setFromTriplets(triplets.begin(), triplets.end());

    return globalStiffnessMatrix;
}

Eigen::MatrixXd Mesh::parallelAssembleStiffnessMatrix() {
    Eigen::MatrixXd globalStiffnessMatrix(Nodes.size() * 2, Nodes.size() * 2);
    globalStiffnessMatrix.setZero();

    std::mutex mtx;
    auto task = [&](Element* element) {
        std::lock_guard<std::mutex> lock(mtx);
        auto&& ke = element->getStiffnessMatrix();
        for (size_t i = 0; i < Serendipity::nodeNum; ++i) {
            for (size_t j = 0; j < Serendipity::nodeNum; ++j) {
                globalStiffnessMatrix(2 * element->nodes[i]->getIndex(),
                                      2 * element->nodes[j]->getIndex()) +=
                    ke(2 * i, 2 * j);
                globalStiffnessMatrix(2 * element->nodes[i]->getIndex(),
                                      2 * element->nodes[j]->getIndex() + 1) +=
                    ke(2 * i, 2 * j + 1);
                globalStiffnessMatrix(2 * element->nodes[i]->getIndex() + 1,
                                      2 * element->nodes[j]->getIndex()) +=
                    ke(2 * i + 1, 2 * j);
                globalStiffnessMatrix(2 * element->nodes[i]->getIndex() + 1,
                                      2 * element->nodes[j]->getIndex() + 1) +=
                    ke(2 * i + 1, 2 * j + 1);
            }
        }
    };

    ThreadPool pool(std::thread::hardware_concurrency());
    for (auto element : Elements) {
        pool.enqueue([&task, &element] { return task(element); });
    }

    return globalStiffnessMatrix;
}
/*
Eigen::SparseMatrix<double> Mesh::parallelSparseAssembleStiffnessMatrix() {
    Eigen::SparseMatrix<double> globalStiffnessMatrix(Nodes.size() * 2,
                                                      Nodes.size() * 2);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(Elements.size() * 4 * Serendipity::nodeNum *
                     Serendipity::nodeNum);

    std::mutex mtx;
    auto threadTask = [&](Element* element) {
        auto&& ke = element->stiffnessMatrix();
        for (size_t i = 0; i < Serendipity::nodeNum; ++i) {
            for (size_t j = 0; j < Serendipity::nodeNum; ++j) {
                std::lock_guard<std::mutex> lock(mtx);
                triplets.push_back({2 * element->nodes[i]->getIndex(),
                                    2 * element->nodes[j]->getIndex(),
                                    ke(2 * i, 2 * j)});
                triplets.push_back({2 * element->nodes[i]->getIndex(),
                                    2 * element->nodes[j]->getIndex() + 1,
                                    ke(2 * i, 2 * j + 1)});
                triplets.push_back({2 * element->nodes[i]->getIndex() + 1,
                                    2 * element->nodes[j]->getIndex(),
                                    ke(2 * i + 1, 2 * j)});
                triplets.push_back({2 * element->nodes[i]->getIndex() + 1,
                                    2 * element->nodes[j]->getIndex() + 1,
                                    ke(2 * i + 1, 2 * j + 1)});
            }
        }
    };

    ThreadPool pool(std::thread::hardware_concurrency());
    for (auto element : Elements) {
        pool.enqueue([&threadTask, &element] { return threadTask(element);
        });
    }

    globalStiffnessMatrix.setFromTriplets(triplets.begin(), triplets.end());

    return globalStiffnessMatrix;
}
 */
int Mesh::Solve(std::list<Load>& loads, std::list<Boundary>& boundaries,
                bool verbose) {
    Force.resize(Nodes.size() * 2);
    Force.setZero();
    for (auto&& load : loads) {
        auto&& equivalentForce = Mesh::equivalentForce(&load);
        for (size_t i = 0; i < 3; ++i) {
            Force.coeffRef(2 * load.nodes[i]->getIndex()) +=
                equivalentForce[2 * i];
            Force.coeffRef(2 * load.nodes[i]->getIndex() + 1) +=
                equivalentForce[2 * i + 1];
        }
    }

    Timer timer;

    auto&& K = sparseAssembleStiffnessMatrix();

    if (verbose) {
        std::cout << "  Stiffness matrix assembled in " << timer.elapsed()
                  << " ms" << std::endl;

        // Apply boundary conditions
        timer.reset();
    }

    Eigen::VectorXd List1 = Eigen::VectorXd::Ones(Nodes.size() * 2);
    Eigen::VectorXd List2 = Eigen::VectorXd::Zero(Nodes.size() * 2);
    for (auto&& boundary : boundaries) {
        for (size_t i = 0; i < 2; ++i) {
            if (boundary.fixed[i]) {
                size_t j = 2 * boundary.node->getIndex() + i;
                List1(j) = 0;
                List2(j) = 1;
                Force.coeffRef(j) = 0;
            }
        }
    }
    {
        Eigen::SparseMatrix<double> diagonalMatrix1;
        Eigen::SparseMatrix<double> diagonalMatrix2;
        diagonalMatrix1 = List1.asDiagonal();
        diagonalMatrix2 = List2.asDiagonal();
        K = diagonalMatrix1 * K * diagonalMatrix1 + diagonalMatrix2;
    }
    /*
        for (auto&& boundary : boundaries) {
            for (size_t i = 0; i < 2; ++i) {
                if (boundary.fixed[i]) {
                    size_t j = 2 * boundary.node->getIndex() + i;
                    globalStiffnessMatrix(j, j) *= 1e50;
                    Force(j) = 0;
                }
            }
        }
     */

    if (verbose) {
        std::cout << "  Boundary conditions applied in " << timer.elapsed()
                  << " ms" << std::endl;
    }

    // auto&& K = globalStiffnessMatrix.sparseView();
    Eigen::SparseVector<double> U;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

    if (verbose) {
        timer.reset();
    }

    solver.analyzePattern(K);
    solver.factorize(K);

    if (verbose) {
        std::cout << "  Analyzed pattern and factorized in " << timer.elapsed()
                  << " ms" << std::endl;
    }

    if (solver.info() != Eigen::Success) {
        std::cerr << "LU decomposition failed" << std::endl;
        std::cerr << "Error code: " << solver.info() << std::endl;
        std::cerr << "Error message: " << solver.lastErrorMessage() << "\n";
        return -1;
    }

    if (verbose) {
        timer.reset();
    }

    U = solver.solve(Force);
    if (verbose) {
        std::cout << "  Solver solved in " << timer.elapsed() << " ms"
                  << std::endl;
    }

    for (size_t i = 0; i < Nodes.size(); ++i) {
        Nodes[i]->Displacement(0) = U.coeff(2 * i);
        Nodes[i]->Displacement(1) = U.coeff(2 * i + 1);
    }

    for (auto&& element : Elements) {
        element->calculateStrainStressGaussPoint();
    }

    return 0;
}