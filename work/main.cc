#include <iostream>
#include <toml.hpp>
#include <vector>

#include "ApplyBoundary.h"
#include "GenerateMesh.h"
#include "Material.h"
#include "Mesh.h"
#include "SetMaterial.h"
#include "vtkManager.h"

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " [path to your settings.toml]"
                  << std::endl;
        return 0;
    }
    try {
        auto settings = toml::parse_file(argv[1]);
        // Define geometry
        double L = settings["Rectangle"]["L"].value_or(2);
        double B = settings["Rectangle"]["B"].value_or(0.9);
        double ksi = settings["Ellipse"]["ksi"].value_or(3);
        double a = B / ksi;
        double b = settings["Ellipse"]["b"].value_or(0.1);
        double lc = settings["Mesh"]["size"].value_or(0.02);
        double rf = settings["Mesh"]["refinementFactor"].value_or(8);
        bool isSerendipity = settings["Mesh"]["isSerendipity"].value_or(true);
        int Algorithm = settings["Mesh"]["Algorithm"].value_or(8);

        // Generate mesh
        std::vector<double> nodeCoord;
        std::vector<size_t> elementNodeTags;
        int err = generate_mesh(nodeCoord, elementNodeTags, L, B, a, b, lc, rf,
                                isSerendipity, Algorithm);
        if (err != 0) {
            return err;
        }
        std::cout << "Creating mesh..." << std::endl;
        Mesh mesh(Mesh::MeshType::serendipity, nodeCoord, elementNodeTags);
        std::cout << "Mesh created with " << mesh.Nodes.size() << " nodes and "
                  << mesh.Elements.size() << " elements" << std::endl;
        nodeCoord.clear();
        elementNodeTags.clear();

        // Set material
        double E1 = settings["Material"]["Base"][0].value_or(1);
        double nu1 = settings["Material"]["Base"][1].value_or(0.3);
        double E2 = settings["Material"]["Inclusion"][1].value_or(1);
        double nu2 = settings["Material"]["Inclusion"][1].value_or(0.2);
        Elastic base(E1, nu1);
        Elastic inclusion(E2, nu2);
        set_material(&mesh, {&base, &inclusion}, a, b);

        // Apply boundary
        double value = settings["Load"]["Value"].value_or(1);
        auto&& loadCondition = apply_load(&mesh, L, B, value);
        auto&& boundaryCondition = apply_boundary(&mesh, 0, 0);

        // for (auto&& element : mesh.Elements) {
        //     std::cout << "Element " << element->index
        //               << " stiffness matrix:" << std::endl;
        //     std::cout << element->stiffnessMatrix() << std::endl;
        // }

        mesh.Solve(loadCondition, boundaryCondition);

        vtkManager vtk(mesh);
        vtk.setData(mesh);
        vtk.write("result.vtk");
    } catch (const std::exception& e) {
        std::cerr << e.what();
        return -1;
    }

    return 0;
}