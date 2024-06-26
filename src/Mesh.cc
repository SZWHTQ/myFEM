#include <omp.h>

#include <Eigen/SparseLU>
#include <climits>
#include <iostream>
#include <memory>
#include <mutex>
#include <new>
#include <thread>
#include <vector>

#include "Boundary.h"
#include "Element.h"
#include "Mesh.h"
#include "Serendipity.h"
#include "ThreadPool.h"
#include "Timer.h"

Mesh::Mesh(const MeshType type, const std::vector<double> nodeCoord,
           const std::vector<size_t> elementNodeTags,
           const std::vector<size_t>& boundaryNodeTags, const bool planeStress_)
    : planeStress(planeStress_) {
    nodeNum = nodeCoord.size() / 3;
    Nodes.reserve(nodeNum);
    for (size_t i = 0; i < nodeCoord.size(); i += 3) {
        size_t id = i / 3;
        Nodes.push_back(std::make_shared<Node>(
            id, nodeCoord[i], nodeCoord[i + 1], nodeCoord[i + 2]));
    }

    for (unsigned long boundaryNodeTag : boundaryNodeTags) {
        Nodes[boundaryNodeTag]->isBoundary = true;
    }

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

Mesh::Mesh(const std::vector<double> nodeCoord,
           const std::vector<std::pair<MeshType, std::vector<size_t>>>
               elementsTypeAndNodeTags,
           const std::vector<size_t>& boundaryNodeTags, const bool planeStress_)
    : planeStress(planeStress_) {
    nodeNum = nodeCoord.size() / 3;
    Nodes.reserve(nodeNum);
    for (size_t i = 0; i < nodeCoord.size(); i += 3) {
        Nodes.push_back(std::make_shared<Node>(
            i / 3, nodeCoord[i], nodeCoord[i + 1], nodeCoord[i + 2]));
    }

    for (unsigned long boundaryNodeTag : boundaryNodeTags) {
        Nodes[boundaryNodeTag - 1]->isBoundary = true;
    }

    for (auto&& [type, elementNodeTags] : elementsTypeAndNodeTags) {
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

Mesh::~Mesh() {
    for (Element* elem : Elements) {
        delete elem;
    }
    Elements.clear();
}

std::vector<double> const Mesh::equivalentForce(const Load& load) {
    std::vector<double> equivalentForce(6);
    Eigen::VectorXd X(6), Y(6);
    auto& n = load.nodes;
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
        (X(0) * load.shearForce[0] + X(1) * load.shearForce[1] +
         X(2) * load.shearForce[2] + Y(0) * load.normalForce[0] +
         Y(1) * load.normalForce[1] + Y(2) * load.normalForce[2]) *
        load.thickness / 30;
    equivalentForce[1] =
        (Y(0) * load.shearForce[0] + Y(1) * load.shearForce[1] +
         Y(2) * load.shearForce[2] - X(0) * load.normalForce[0] -
         X(1) * load.normalForce[1] - X(2) * load.normalForce[2]) *
        load.thickness / 30;
    equivalentForce[2] =
        (X(1) * load.shearForce[0] + X(3) * load.shearForce[1] +
         X(4) * load.shearForce[2] + Y(1) * load.normalForce[0] +
         Y(3) * load.normalForce[1] + Y(4) * load.normalForce[2]) *
        load.thickness / 30;
    equivalentForce[3] =
        (Y(1) * load.shearForce[0] + Y(3) * load.shearForce[1] +
         Y(4) * load.shearForce[2] - X(1) * load.normalForce[0] -
         X(3) * load.normalForce[1] - X(4) * load.normalForce[2]) *
        load.thickness / 30;
    equivalentForce[4] =
        (X(2) * load.shearForce[0] + X(4) * load.shearForce[1] +
         X(5) * load.shearForce[2] + Y(2) * load.normalForce[0] +
         Y(4) * load.normalForce[1] + Y(5) * load.normalForce[2]) *
        load.thickness / 30;
    equivalentForce[5] =
        (Y(2) * load.shearForce[0] + Y(4) * load.shearForce[1] +
         Y(5) * load.shearForce[2] - X(2) * load.normalForce[0] -
         X(4) * load.normalForce[1] - X(5) * load.normalForce[2]) *
        load.thickness / 30;
    return equivalentForce;
}

Eigen::MatrixXd Mesh::assembleStiffnessMatrix() {
    Eigen::MatrixXd globalStiffnessMatrix(Nodes.size() * 2, Nodes.size() * 2);
    globalStiffnessMatrix.setZero();

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

Eigen::SparseMatrix<double> Mesh::assembleSparseStiffnessMatrix() {
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
    {
        ThreadPool pool(std::thread::hardware_concurrency());
        for (const auto element : Elements) {
            pool.enqueue([&task, &element] { return task(element); });
        }
    }

    return globalStiffnessMatrix;
}

// FIXME: This function is wrong
Eigen::SparseMatrix<double> Mesh::parallelAssembleSparseStiffnessMatrix() {
    Eigen::SparseMatrix<double> globalStiffnessMatrix(Nodes.size() * 2,
                                                      Nodes.size() * 2);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(Elements.size() * 4 * Serendipity::nodeNum *
                     Serendipity::nodeNum);

    // std::mutex mtx;
    auto threadTask = [&](Element* element) {
        auto&& ke = element->getStiffnessMatrix();
        for (size_t i = 0; i < Serendipity::nodeNum; ++i) {
            for (size_t j = 0; j < Serendipity::nodeNum; ++j) {
                // std::lock_guard<std::mutex> lock(mtx);
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
    };
    {
        ThreadPool pool(std::thread::hardware_concurrency());
        for (const auto element : Elements) {
            pool.enqueue(
                [&threadTask, &element] { return threadTask(element); });
        }
    }

    globalStiffnessMatrix.setFromTriplets(triplets.begin(), triplets.end());

    return globalStiffnessMatrix;
}

int Mesh::Solve(const std::list<Load>& loads,
                const std::list<Boundary>& boundaries, const bool verbose) {
    Timer timer;
    Eigen::setNbThreads(omp_get_max_threads());
    if (verbose) {
        unsigned int eigenThreadsNum = Eigen::nbThreads();
        std::cout << "  Eigen threads: " << eigenThreadsNum << std::endl;
    }

    // Get equivalent force
    this->Force.resize(Nodes.size() * 2);
    this->Force.setZero();
    for (const auto& load : loads) {
        auto& equivalentForce = Mesh::equivalentForce(load);
        for (size_t i = 0; i < 3; ++i) {
            this->Force.coeffRef(2 * load.nodes[i]->getIndex()) +=
                equivalentForce[2 * i];
            this->Force.coeffRef(2 * load.nodes[i]->getIndex() + 1) +=
                equivalentForce[2 * i + 1];
        }
    }
    if (verbose) {
        std::cout << "  Equivalent force calculated in " << timer << std::endl;
        timer.reset();
    }

    // Assemble stiffness matrix
    auto K = assembleSparseStiffnessMatrix();
    if (verbose) {
        std::cout << "  Stiffness matrix assembled in " << timer << std::endl;
        timer.reset();
    }

    // Apply boundary conditions
    {
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
        std::cout << "  Boundary conditions applied in " << timer << std::endl;
    }

    // Solve
    // auto&& K = globalStiffnessMatrix.sparseView();
    Eigen::SparseVector<double> U;
    try {
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        if (verbose) {
            timer.reset();
        }
        solver.analyzePattern(K);
        solver.factorize(K);
        if (verbose) {
            std::cout << "  Analyzed pattern and factorized in " << timer
                      << std::endl;
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
    } catch (const std::bad_alloc& e) {
        std::cerr << "Memory allocation failed" << std::endl;
        std::cerr << "Use Conjugate Gradient instead" << std::endl;

        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                                 Eigen::Lower | Eigen::Upper,
                                 Eigen::DiagonalPreconditioner<double>>
            solver;
        if (verbose) {
            timer.reset();
        }

        solver.compute(K);

        solver.setTolerance(1e-9);
        solver.setMaxIterations(INT_MAX);
        U = solver.solve(Force);

        if (solver.info() != Eigen::Success) {
            std::cerr << "Conjugate Gradient failed" << std::endl;
            std::cerr << "Error code: " << solver.info() << std::endl;
            return -1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Solver failed" << std::endl;
        std::cerr << e.what() << std::endl;
        return -1;
    }
    if (verbose) {
        std::cout << "  Solver solved in " << timer << std::endl;
    }

    // Copy displacement to Nodes
    for (size_t i = 0; i < Nodes.size(); ++i) {
        Nodes[i]->Displacement(0) = U.coeff(2 * i);
        Nodes[i]->Displacement(1) = U.coeff(2 * i + 1);
    }

    // Calculate stain stress
    for (auto&& element : Elements) {
        element->calculateNodeStrainStressViaGaussPoint();
    }

    return 0;
}
