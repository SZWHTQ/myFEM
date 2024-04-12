#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkQuadraticQuad.h>

#include "Material.h"
#include "vtkManager.h"

vtkManager::vtkManager(Mesh& mesh) {
    vtkNew<vtkPoints> Points;
    for (auto&& node : mesh.Nodes) {
        Points->InsertNextPoint(node->x, node->y, node->z);
    }
    Grid->SetPoints(Points);
    for (auto&& element : mesh.Elements) {
        vtkNew<vtkQuadraticQuad> QuadraticQuad;
        for (size_t i = 0; i < element->nodes.size(); ++i) {
            QuadraticQuad->GetPointIds()->SetId(i,
                                                element->nodes[i]->getIndex());
        }
        Grid->InsertNextCell(QuadraticQuad->GetCellType(),
                             QuadraticQuad->GetPointIds());
    }
};

void vtkManager::setMeshData(Mesh& mesh) {
    vtkNew<vtkFloatArray> Displacement;
    Displacement->SetNumberOfComponents(3);
    Displacement->SetName("Displacement");

    for (auto&& node : mesh.Nodes) {
        Displacement->InsertNextTuple(node->Displacement.data());
    }
    Grid->GetPointData()->AddArray(Displacement);

    vtkNew<vtkFloatArray> Force;
    Force->SetNumberOfComponents(3);
    Force->SetName("Force");
    for (size_t i = 0; i < mesh.Force.size() / 2; ++i) {
        Force->InsertNextTuple3(mesh.Force.coeff(i * 2),
                                mesh.Force.coeff(i * 2 + 1), 0);
    }
    Grid->GetPointData()->AddArray(Force);

    vtkNew<vtkFloatArray> Normal;
    Normal->SetNumberOfComponents(3);
    Normal->SetName("FaceNormal");
    for (auto&& element : mesh.Elements) {
        Eigen::Vector3d edge1{element->nodes[1]->x - element->nodes[0]->x,
                              element->nodes[1]->y - element->nodes[0]->y,
                              element->nodes[1]->z - element->nodes[0]->z};
        Eigen::Vector3d edge2{element->nodes[2]->x - element->nodes[0]->x,
                              element->nodes[2]->y - element->nodes[0]->y,
                              element->nodes[2]->z - element->nodes[0]->z};
        Eigen::Vector3d normal = edge1.cross(edge2);
        normal.normalize();
        Normal->InsertNextTuple(normal.data());
    }
    Grid->GetCellData()->AddArray(Normal);

    vtkNew<vtkFloatArray> Strain;
    Strain->SetNumberOfComponents(3);
    Strain->SetName("Strain");
    for (auto&& node : mesh.Nodes) {
        Strain->InsertNextTuple(node->Strain.data());
    }
    Grid->GetPointData()->AddArray(Strain);

    vtkNew<vtkFloatArray> Stress;
    Stress->SetNumberOfComponents(3);
    Stress->SetName("Stress");
    for (auto&& node : mesh.Nodes) {
        Stress->InsertNextTuple(node->Stress.data());
    }
    Grid->GetPointData()->AddArray(Stress);

    vtkNew<vtkFloatArray> Material;
    Material->SetNumberOfComponents(1);
    Material->SetName("Material");
    for (auto&& element : mesh.Elements) {
        Material->InsertNextValue(element->material->getIndex());
    }
    Grid->GetCellData()->AddArray(Material);

    vtkNew<vtkFloatArray> AreaArray;
    AreaArray->SetNumberOfComponents(1);
    AreaArray->SetName("Area");
    for (auto&& element : mesh.Elements) {
        AreaArray->InsertNextValue(element->getArea());
    }
    Grid->GetCellData()->AddArray(AreaArray);

    vtkNew<vtkFloatArray> strainEnergy;
    strainEnergy->SetNumberOfComponents(1);
    strainEnergy->SetName("StrainEnergy");
    for (auto&& element : mesh.Elements) {
        strainEnergy->InsertNextValue(element->getStrainEnergy());
    }
    Grid->GetCellData()->AddArray(strainEnergy);

    vtkNew<vtkIntArray> Boundary;
    Boundary->SetNumberOfComponents(1);
    Boundary->SetName("Boundary");
    for (auto&& node : mesh.Nodes) {
        Boundary->InsertNextValue(node->isBoundary);
    }
    Grid->GetPointData()->AddArray(Boundary);

}

void vtkManager::write(std::string fileName, bool isBinary) {
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(Grid);
    if (isBinary) {
        writer->SetDataModeToBinary();
    } else {
        writer->SetDataModeToAscii();
    }
    writer->Write();
}