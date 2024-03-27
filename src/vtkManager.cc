#include <vtkFloatArray.h>
#include <vtkManager.h>
#include <vtkNew.h>
#include <vtkPointData.h>

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

void vtkManager::setData(std::string name, Mesh& mesh) {
    vtkNew<vtkFloatArray> displacement;
    displacement->SetNumberOfComponents(3);
    displacement->SetName(name.c_str());

    for (auto&& node : mesh.Nodes) {
        displacement->InsertNextTuple3(
            node.displacement(0), node.displacement(1), node.displacement(2));
    }
    Grid->GetPointData()->SetVectors(displacement);
}

void vtkManager::write(std::string fileName) {
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(Grid);
    writer->Write();
}