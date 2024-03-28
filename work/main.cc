#include <iostream>
#include <vector>

#include "ApplyBoundary.h"
#include "GenerateMesh.h"
#include "Material.h"
#include "Mesh.h"
#include "SetMaterial.h"
#include "vtkManager.h"

int main() {
    // Define geometry
    double L = 2;
    double B = 0.9;
    double ksi = 3;
    double a = 0.9 / ksi;
    double b = 0.1;
    double lc = 0.02;

    // Generate mesh
    std::vector<double> nodeCoord;
    std::vector<size_t> elementNodeTags;
    int err = generate_mesh(nodeCoord, elementNodeTags, L, B, a, b, lc);
    if (err != 0) {
        return -1;
    }
    Mesh mesh(Mesh::MeshType::serendipity, nodeCoord, elementNodeTags);
    std::cout << "Mesh created with " << mesh.Nodes.size() << " nodes and "
              << mesh.Elements.size() << " elements" << std::endl;

    // Set material
    Elastic base(1, 0.3);
    Elastic inclusion(1, 0.3);
    set_material(&mesh, {&base, &inclusion}, a, b);

    // Apply boundary
    auto&& loadCondition = apply_load(&mesh, L, B);
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

    return 0;
}