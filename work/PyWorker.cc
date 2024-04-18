#if _WIN32 || _WIN64
#define EXPORT __declspec(dllexport)
#else
#define EXPORT __attribute__((visibility("default")))
#endif

#include <gmsh.h>

#include <iomanip>
#include <iostream>
#include <vector>

#include "ApplyBoundary.h"
#include "GenerateMesh.h"
#include "GetStrainEnergyChange.h"
#include "Material.h"
#include "Mesh.h"
#include "SetMaterial.h"
#include "Timer.h"

#ifdef __cplusplus
extern "C" {
#endif
EXPORT double PyWorker(const double L, const double B, const double ksi,
                       const double a_b, const double loadValue,
                       const double matrixElasticModulus,
                       const double matrixPoissonRatio,
                       const double inclusionElasticModulus,
                       const double inclusionPoissonRatio,
                       const double meshSize, const double refinementFactor,
                       const int meshAlgorithm, const bool isSerendipity,
                       const bool convertToSquare, const bool writeInp,
                       const bool writeMsh, const bool runFltk,
                       const bool isPlaneStress, const bool verbose) {
    try {
        // Define geometry
        double a = B / ksi;
        double b = a / a_b;
        if (verbose) {
            std::cout << "L: " << L << " B: " << B << std::endl;
            std::cout << "a: " << a << " b: " << b << std::endl;
            std::cout << "matrix elastic modulus: " << matrixElasticModulus
                      << " matrix poisson ratio: " << matrixPoissonRatio
                      << std::endl;
            std::cout << "inclusion elastic modulus: "
                      << inclusionElasticModulus
                      << " inclusion poisson ratio: " << inclusionPoissonRatio
                      << std::endl;
        }

        // Generate mesh
        Timer t, timer;
        std::vector<double> nodeCoord;
        std::vector<size_t> elementNodeTags;
        std::vector<size_t> elementMaterialTags;
        std::vector<size_t> interfaceNodeTags;

        // Initialize the Gmsh library
        gmsh::initialize();
        int err = generate_mesh(
            nodeCoord, elementNodeTags, elementMaterialTags, interfaceNodeTags,
            L, B, a, b, meshSize, refinementFactor, isSerendipity,
            meshAlgorithm, convertToSquare, writeInp, writeMsh, runFltk);
        if (err != 0) {
            return err;
        }
        if (verbose) {
            std::cout << std::fixed << std::setprecision(2);
            std::cout << "Mesh created in " << timer.elapsed() << " ms"
                      << std::endl;
        }

        // Convert mesh
        if (verbose) {
            timer.reset();
        }
        Mesh mesh(nodeCoord,
                  {std::pair(Mesh::MeshType::serendipity, elementNodeTags)},
                  interfaceNodeTags, isPlaneStress);
        if (verbose) {
            std::cout << "Mesh converted in " << timer.elapsed() << " ms"
                      << " with " << mesh.Nodes.size() << " nodes and "
                      << mesh.Elements.size() << " elements" << std::endl;
        }

        // Set material
        Elastic matrix(1, matrixElasticModulus, matrixPoissonRatio);
        Elastic inclusion(2, inclusionElasticModulus, inclusionPoissonRatio);
        // set_material(&mesh, {&matrix, &inclusion}, a, b);
        set_material(&mesh, {&matrix, &inclusion}, elementMaterialTags);

        // Apply boundary
        auto&& loadCondition = apply_load(&mesh, L, B, loadValue);
        auto&& boundaryCondition = apply_boundary(&mesh, 0, 0);

        // Solve
        if (verbose) {
            timer.reset();
        }
        mesh.Solve(loadCondition, boundaryCondition, false);
        if (verbose) {
            std::cout << "Mesh solved in " << timer.elapsed() << " ms"
                      << std::endl;
        }

        double matrixStrainEnergy = 0, inclusionStrainEnergy = 0;
        for (auto& element : mesh.Elements) {
            if (element->material->getIndex() == 1) {
                matrixStrainEnergy += element->getStrainEnergy();
            } else if (element->material->getIndex() == 2) {
                inclusionStrainEnergy += element->getStrainEnergy();
            } else {
                std::cerr << "Wrong material index: "
                          << element->material->getIndex() << " at "
                          << element->getIndex() << std::endl;
            }
        }

        auto deltaU =
            getStrainEnergyChange(&mesh, &matrix, &inclusion, isPlaneStress);
        std::cout << std::setprecision(-1);
        if (verbose) {
            std::cout << "Strain energy change: " << deltaU * 4 << std::endl;
        }
        return deltaU * 4;

        // Clear memory
        {
            nodeCoord.clear();
            std::vector<double>().swap(nodeCoord);
            elementNodeTags.clear();
            std::vector<size_t>().swap(elementNodeTags);
        }

        // Finalize the Gmsh library
        gmsh::finalize();

    } catch (const std::exception& e) {
        std::cerr << e.what();
        return -1;
    }
}
#ifdef __cplusplus
}
#endif