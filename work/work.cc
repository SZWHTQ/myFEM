#include <gmsh.h>

#include <iostream>
#include <vector>

#include "ApplyBoundary.h"
#include "GenerateMesh.h"
#include "GetStrainEnergyChange.h"
#include "Material.h"
#include "Mesh.h"
#include "SetMaterial.h"
#include "Timer.h"
#include "toml.hpp"
#include "work.h"

extern "C" {
EXPORT double work(const char* tomlFilePath,
                   const double inclusionElasticModulus, const double ksi,
                   const double a_b, const bool verbose) {
    std::cout << inclusionElasticModulus << std::endl;
    std::cout << ksi << std::endl;
    std::cout << a_b << std::endl;
    try {
        auto settings = toml::parse_file(tomlFilePath);
        // Define geometry
        double L = settings["Rectangle"]["L"].value_or(2);
        double B = settings["Rectangle"]["B"].value_or(0.9);
        double a = B / ksi;
        double b = a / a_b;
        bool convertToSquare =
            settings["Ellipse"]["convertToSquare"].value_or(false);
        double lc = settings["Mesh"]["size"].value_or(0.02);
        double refinementFactor =
            settings["Mesh"]["refinementFactor"].value_or(8);
        int meshAlgorithm = settings["Mesh"]["Algorithm"].value_or(8);
        bool isSerendipity = settings["Mesh"]["Serendipity"].value_or(true);
        bool isPlaneStress = settings["Mesh"]["planeStress"].value_or(true);

        // Generate mesh
        Timer t, timer;
        std::vector<double> nodeCoord;
        std::vector<size_t> elementNodeTags;
        std::vector<size_t> elementMaterialTags;
        std::vector<size_t> interfaceNodeTags;

        // Initialize the Gmsh library
        gmsh::initialize();
        int err =
            generate_mesh(nodeCoord, elementNodeTags, elementMaterialTags,
                          interfaceNodeTags, L, B, a, b, lc, refinementFactor,
                          isSerendipity, meshAlgorithm, convertToSquare);
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
        double E1 = settings["Material"]["Matrix"][0].value_or(1.0);
        double nu1 = settings["Material"]["Matrix"][1].value_or(0.3);
        double E2 = inclusionElasticModulus;
        double nu2 = settings["Material"]["Inclusion"][1].value_or(0.2);
        Elastic matrix(1, E1, nu1);
        Elastic inclusion(2, E2, nu2);
        // set_material(&mesh, {&matrix, &inclusion}, a, b);
        set_material(&mesh, {&matrix, &inclusion}, elementMaterialTags);

        // Apply boundary
        double value = settings["Load"]["Value"].value_or(1);
        auto&& loadCondition = apply_load(&mesh, L, B, value);
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
}