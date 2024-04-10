#include <gmsh.h>

#include <iostream>
#include <toml.hpp>
#include <vector>

#include "ApplyBoundary.h"
#include "GenerateMesh.h"
#include "GetStrainEnergyChange.h"
#include "Material.h"
#include "Mesh.h"
#include "SetMaterial.h"
#include "Timer.h"
#include "vtkManager.h"

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " [path to your settings.toml]"
                  << std::endl;
        return 0;
    }
    try {
        std::cout << "Parsing..." << std::endl;

        auto settings = toml::parse_file(argv[1]);
        // Define geometry
        double L = settings["Rectangle"]["L"].value_or(2);
        double B = settings["Rectangle"]["B"].value_or(0.9);
        double ksi = settings["Ellipse"]["ksi"].value_or(3);
        double a = B / ksi;
        double a_b = settings["Ellipse"]["a_b"].value_or(1. / 3);
        double b = a / a_b;
        double lc = settings["Mesh"]["size"].value_or(0.02);
        double rf = settings["Mesh"]["refinementFactor"].value_or(8);
        bool isSerendipity = settings["Mesh"]["Serendipity"].value_or(true);
        int Algorithm = settings["Mesh"]["Algorithm"].value_or(8);
        bool isPlaneStress = settings["Mesh"]["planeStress"].value_or(true);

        // Generate mesh
        Timer t, timer;
        std::vector<double> nodeCoord;
        std::vector<size_t> elementNodeTags;
        std::vector<size_t> boundaryNodeTags;

        int err = generate_mesh(nodeCoord, elementNodeTags, boundaryNodeTags, L,
                                B, a, b, lc, rf, isSerendipity, Algorithm);
        if (err != 0) {
            return err;
        }
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "Mesh created in " << timer.elapsed() << " ms"
                  << std::endl;

        // Convert mesh
        timer.reset();
        Mesh mesh(nodeCoord,
                  {std::pair(Mesh::MeshType::serendipity, elementNodeTags)},
                  boundaryNodeTags, isPlaneStress);
        std::cout << "Mesh converted in " << timer.elapsed() << " ms"
                  << " with " << mesh.Nodes.size() << " nodes and "
                  << mesh.Elements.size() << " elements" << std::endl;

        // Set material
        double E1 = settings["Material"]["Matrix"][0].value_or(1);
        double nu1 = settings["Material"]["Matrix"][1].value_or(0.3);
        double E2 = settings["Material"]["Inclusion"][0].value_or(1);
        double nu2 = settings["Material"]["Inclusion"][1].value_or(0.2);
        Elastic matrix(1, E1, nu1);
        Elastic inclusion(2, E2, nu2);
        set_material(&mesh, {&matrix, &inclusion}, a, b);

        // Apply boundary
        double value = settings["Load"]["Value"].value_or(1);
        auto&& loadCondition = apply_load(&mesh, L, B, value);
        auto&& boundaryCondition = apply_boundary(&mesh, 0, 0);

        // Solve
        timer.reset();
        mesh.Solve(loadCondition, boundaryCondition, true);
        std::cout << "Mesh solved in " << timer.elapsed() << " ms" << std::endl;

        // Write output
        std::string vtkFileName =
            settings["VTK"]["fileName"].value_or("result.vtk");
        bool isBinary = settings["VTK"]["Binary"].value_or(false);
        timer.reset();
        vtkManager vtk(mesh);
        vtk.setMeshData(mesh);
        vtk.write(vtkFileName, isBinary);
        std::cout << "Vtk written in " << timer.elapsed() << " ms" << std::endl;
        std::cout << "Total time: " << t.elapsed() << " ms" << std::endl;
        std::cout << "\n";

        double matrixStrainEnergy = 0, inclusionStrainEnergy = 0;
        for (auto& element : mesh.Elements) {
            if (element->material->getIndex() == 1) {
                matrixStrainEnergy += element->getStrainEnergy();
            } else if (element->material->getIndex() == 2) {
                inclusionStrainEnergy += element->getStrainEnergy();
            } else {
                std::cerr << "Wrong material index" << std::endl;
            }
        }

        std::cout << std::fixed << std::setprecision(-1);
        std::cout << "Matrix with inclusion:" << std::endl;
        std::cout << "  Matrix strain energy: " << matrixStrainEnergy
                  << std::endl;
        std::cout << "  Inclusion strain energy: " << inclusionStrainEnergy
                  << std::endl;
        std::cout << std::fixed << std::setprecision(2);

        Mesh meshNoInclusion(
            nodeCoord,
            {std::pair(Mesh::MeshType::serendipity, elementNodeTags)},
            boundaryNodeTags, isPlaneStress);
        Material temp(2, matrix.getElasticModulus(), matrix.getPoissonRatio());
        set_material(&meshNoInclusion, {&matrix, &temp}, a, b);
        meshNoInclusion.Solve(loadCondition, boundaryCondition);
        vtkManager vtkNoInclusion(meshNoInclusion);
        vtkNoInclusion.setMeshData(meshNoInclusion);
        vtkNoInclusion.write("PlateNoInclusion.vtu", isBinary);

        double matrixStrainEnergyNoInclusion = 0,
               inclusionStrainEnergyNoInclusion = 0;
        for (auto& element : meshNoInclusion.Elements) {
            if (element->material->getIndex() == 1) {
                matrixStrainEnergyNoInclusion += element->getStrainEnergy();
            } else if (element->material->getIndex() == 2) {
                inclusionStrainEnergyNoInclusion += element->getStrainEnergy();
            } else {
                std::cerr << "Wrong material index" << std::endl;
            }
        }

        std::cout << std::fixed << std::setprecision(-1);
        std::cout << "Matrix no inclusion:" << std::endl;
        std::cout << "  Matrix strain energy: " << matrixStrainEnergyNoInclusion
                  << std::endl;
        std::cout << "  Inclusion strain energy: "
                  << inclusionStrainEnergyNoInclusion << std::endl;

        // Strain Energy change
        std::cout << "Matrix Strain energy change: "
                  << matrixStrainEnergyNoInclusion - matrixStrainEnergy
                  << std::endl;
        std::cout << "Inclusion Strain energy change: "
                  << inclusionStrainEnergyNoInclusion - inclusionStrainEnergy
                  << std::endl;
        double deltaEnergy =
            matrixStrainEnergyNoInclusion - matrixStrainEnergy +
            inclusionStrainEnergyNoInclusion - inclusionStrainEnergy;
        std::cout << "Strain energy change: " << deltaEnergy << std::endl;
        auto deltaU = getStrainEnergyChange(&mesh, &meshNoInclusion, &matrix,
                                            &inclusion, isPlaneStress);
        std::cout << "Strain energy change: " << deltaU << std::endl;
        std::cout << std::fixed << std::setprecision(2);

        std::cout << "Relative error: "
                  << std::abs(deltaU - deltaEnergy) / deltaEnergy * 100 << "%"
                  << std::endl;

        // Clear memory
        {
            nodeCoord.clear();
            std::vector<double>().swap(nodeCoord);
            elementNodeTags.clear();
            std::vector<size_t>().swap(elementNodeTags);
        }
    } catch (const std::exception& e) {
        std::cerr << e.what();
        return -1;
    }

    gmsh::finalize();

    return 0;
}