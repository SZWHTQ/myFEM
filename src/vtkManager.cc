#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>

#include "vtkManager.h"

vtkManager::vtkManager(Mesh& mesh) {
    vtkNew<vtkPoints> Points;
    for (auto&& node : mesh.Nodes) {
        Points->InsertNextPoint(node.x, node.y, node.z);
    }
    Grid->SetPoints(Points);
    for (auto&& element : mesh.Elements) {
        vtkNew<vtkQuadraticQuad> QuadraticQuad;
        for (int i = 0; i < element->nodes.size(); ++i) {
            QuadraticQuad->GetPointIds()->SetId(i, element->nodes[i]->index);
        }
        Grid->InsertNextCell(QuadraticQuad->GetCellType(),
                             QuadraticQuad->GetPointIds());
    }
};

void vtkManager::setData(Mesh& mesh) {
    vtkNew<vtkFloatArray> Displacement;
    Displacement->SetNumberOfComponents(3);
    Displacement->SetName("Displacement");

    for (auto&& node : mesh.Nodes) {
        Displacement->InsertNextTuple(node.displacement.data());
    }
    Grid->GetPointData()->AddArray(Displacement);

    vtkNew<vtkFloatArray> Force;
    Force->SetNumberOfComponents(3);
    Force->SetName("Force");
    for (size_t i = 0; i < mesh.Force.size() / 2; ++i) {
        Force->InsertNextTuple3(mesh.Force(i * 2), mesh.Force(i * 2 + 1), 0);
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
}

void vtkManager::write(std::string fileName) {
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(Grid);
    writer->Write();
}