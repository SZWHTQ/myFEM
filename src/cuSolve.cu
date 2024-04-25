#include <cuda_runtime_api.h>
#include <cusolverSp.h>
#include <cusparse.h>

#include <cstddef>
#include <iostream>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "Mesh.h"
#include "Timer.h"

// CUDA API error checking
#define CUDA_CHECK(err)                                                   \
    do {                                                                  \
        cudaError_t err_ = (err);                                         \
        if (err_ != cudaSuccess) {                                        \
            printf("CUDA error %d at %s:%d\n", err_, __FILE__, __LINE__); \
            throw std::runtime_error("CUDA error");                       \
        }                                                                 \
    } while (0)

// cusolver API error checking
#define CUSOLVER_CHECK(err)                                                   \
    do {                                                                      \
        cusolverStatus_t err_ = (err);                                        \
        if (err_ != CUSOLVER_STATUS_SUCCESS) {                                \
            printf("cusolver error %d at %s:%d\n", err_, __FILE__, __LINE__); \
            throw std::runtime_error("cusolver error");                       \
        }                                                                     \
    } while (0)

// cusparse API error checking
#define CUSPARSE_CHECK(err)                                                   \
    do {                                                                      \
        cusparseStatus_t err_ = (err);                                        \
        if (err_ != CUSPARSE_STATUS_SUCCESS) {                                \
            printf("cusparse error %d at %s:%d\n", err_, __FILE__, __LINE__); \
            throw std::runtime_error("cusparse error");                       \
        }                                                                     \
    } while (0)

std::unordered_map<Key, double, KeyHash, KeyEqual>
Mesh::getStiffnessMatrixMap() {
    std::unordered_map<Key, double, KeyHash, KeyEqual> stiffnessMatrixMap;
    for (auto&& element : Elements) {
        auto&& ke = element->getStiffnessMatrix();
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                size_t l = element->nodes[i]->getIndex() * 2;
                size_t m = element->nodes[j]->getIndex() * 2;
                stiffnessMatrixMap[Key(l, m)] += ke(2 * i, 2 * j);
                stiffnessMatrixMap[Key(l + 1, m + 1)] +=
                    ke(2 * i + 1, 2 * j + 1);
            }
        }
    }
    return stiffnessMatrixMap;
}

int Mesh::cuSolver(std::list<Load>& loads, std::list<Boundary>& boundaries,
                   bool verbose) {
    Timer timer;

    int numRows = Nodes.size() * 2;

    // Get equivalent force
    Force.resize(Nodes.size() * 2);
    Force.setZero();
    double* hostForce = (double*)(malloc(sizeof(double) * numRows));
    for (auto&& load : loads) {
        auto&& equivalentForce = Mesh::equivalentForce(&load);
        for (size_t i = 0; i < 3; ++i) {
            size_t j = 2 * load.nodes[i]->getIndex();
            Force.coeffRef(j) += equivalentForce[2 * i];
            Force.coeffRef(j + 1) += equivalentForce[2 * i + 1];
            hostForce[j] += equivalentForce[2 * i];
            hostForce[j + 1] += equivalentForce[2 * i + 1];
        }
    }
    if (verbose) {
        std::cout << "  Equivalent force calculated in " << timer << std::endl;
        timer.reset();
    }

    // Assemble stiffness matrix
    auto&& kMap = getStiffnessMatrixMap();
    int nnz = kMap.size();
    if (verbose) {
        std::cout << "  Stiffness matrix assembled in " << timer << std::endl;
        timer.reset();
    }

    // Apply boundary conditions
    {
        for (auto&& boundary : boundaries) {
            for (size_t i = 0; i < 2; ++i) {
                if (boundary.fixed[i]) {
                    size_t j = 2 * boundary.node->getIndex() + i;
                    // kMap[Key(j,j)] = std::numeric_limits<double>::max();
                    kMap[Key(j, j)] *= 1e50;
                    Force.coeffRef(j) = 0;
                    hostForce[j] = 0;
                }
            }
        }
    }

    std::vector<int> hostCsrRowPtr;
    std::vector<int> hostCsrColInd;
    std::vector<double> hostCsrValues;
    int currentRow = 0;
    hostCsrRowPtr.push_back(0);
    for (auto& elem : kMap) {
        int row, col;
        double val;
        std::tie(row, col) = elem.first;
        val = elem.second;

        // Row pointer update
        while (currentRow <= row) {
            hostCsrRowPtr.push_back(hostCsrColInd.size());
            currentRow++;
        }

        // Fill column index and value arrays
        hostCsrColInd.push_back(col);
        hostCsrValues.push_back(val);
    }
    // Deal with the last row
    while (currentRow < numRows) {
        hostCsrRowPtr.push_back(hostCsrColInd.size());
        currentRow++;
    }

    // cuSolver
    cusolverSpHandle_t cusolverSpHandle;
    csrqrInfo_t info = NULL;
    cusparseMatDescr_t descrA = NULL;
    cudaStream_t stream = NULL;
    CUSOLVER_CHECK(cusolverSpCreate(&cusolverSpHandle));
    CUDA_CHECK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    CUSOLVER_CHECK(cusolverSpSetStream(cusolverSpHandle, stream));

    CUSPARSE_CHECK(cusparseCreateMatDescr(&descrA));

    CUSPARSE_CHECK(cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL));
    CUSPARSE_CHECK(cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ONE));

    CUSOLVER_CHECK(cusolverSpCreateCsrqrInfo(&info));

    // Allocate device memory
    double* deviceForce;
    int* deviceCsrRowPtr;
    int* deviceCsrColInd;
    double* deviceCsrValues;
    double* deviceDisplacement;

    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&deviceCsrValues),
                          sizeof(double) * nnz));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&deviceCsrRowPtr),
                          sizeof(int) * (numRows + 1)));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&deviceCsrColInd),
                          sizeof(int) * nnz));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&deviceForce),
                          sizeof(double) * numRows));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&deviceDisplacement),
                          sizeof(double) * numRows));

    int* singularity;
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&singularity), sizeof(int)));

    CUDA_CHECK(cudaMemcpy(deviceCsrValues, hostCsrValues.data(),
                          sizeof(double) * nnz, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(deviceCsrRowPtr, hostCsrRowPtr.data(),
                            sizeof(int) * (numRows + 1), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(deviceCsrColInd, hostCsrColInd.data(),
                            sizeof(int) * nnz, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(deviceForce, hostForce, sizeof(double) * numRows,
                            cudaMemcpyHostToDevice));

    CUSOLVER_CHECK(cusolverSpDcsrlsvluHost(
        cusolverSpHandle, numRows, nnz, descrA, deviceCsrValues,
        deviceCsrRowPtr, deviceCsrColInd, deviceForce, 1e-10, 0,
        deviceDisplacement, singularity));

    double* hostDisplacement = (double*)(malloc(sizeof(double) * numRows));
    CUDA_CHECK(cudaMemcpy(hostDisplacement, deviceDisplacement, sizeof(double) * numRows,
                          cudaMemcpyDeviceToHost));
    
    cudaFree(deviceCsrValues);
    cudaFree(deviceCsrRowPtr);
    cudaFree(deviceCsrColInd);
    cudaFree(deviceForce);
    cudaFree(deviceDisplacement);
    cudaFree(singularity);
    cusolverSpDestroy(cusolverSpHandle);

    
    // Copy displacement to Nodes
    for (size_t i = 0; i < Nodes.size(); ++i) {
        Nodes[i]->Displacement(0) = hostDisplacement[2 * i];
        Nodes[i]->Displacement(1) = hostDisplacement[2 * i + 1];
    }

    // Calculate stain stress
    for (auto&& element : Elements) {
        element->calculateStrainStressGaussPoint();
    }

    return 0;
}  