#include "vtkManager.h"

#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkQuadraticQuad.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include "Element.h"
#include "Material.h"
#include "Mesh.h"

// Implementation class
class vtkManager::vtkManagerImpl {
   public:
    vtkNew<vtkUnstructuredGrid> Grid;
    vtkNew<vtkXMLUnstructuredGridWriter> writer;
};

vtkManager::vtkManager(Mesh& mesh) : pimpl(std::make_unique<vtkManagerImpl>()) {
    vtkNew<vtkPoints> Points;
    for (auto&& node : mesh.Nodes) {
        Points->InsertNextPoint(node->x, node->y, node->z);
    }
    pimpl->Grid->SetPoints(Points);
    for (auto&& element : mesh.Elements) {
        vtkNew<vtkQuadraticQuad> QuadraticQuad;
        for (size_t i = 0; i < element->nodes.size(); ++i) {
            QuadraticQuad->GetPointIds()->SetId(i,
                                                element->nodes[i]->getIndex());
        }
        pimpl->Grid->InsertNextCell(QuadraticQuad->GetCellType(),
                                    QuadraticQuad->GetPointIds());
    }
};

vtkManager::~vtkManager() {
    // Destructor
    // Must be defined in the source file
}

void vtkManager::setMeshData(Mesh& mesh) const {
    vtkNew<vtkFloatArray> Displacement;
    Displacement->SetNumberOfComponents(3);
    Displacement->SetName("Displacement");

    for (auto&& node : mesh.Nodes) {
        Displacement->InsertNextTuple(node->Displacement.data());
    }
    pimpl->Grid->GetPointData()->AddArray(Displacement);

    vtkNew<vtkFloatArray> Force;
    Force->SetNumberOfComponents(3);
    Force->SetName("Force");
    for (size_t i = 0; i < size_t(mesh.Force.size() / 2); ++i) {
        Force->InsertNextTuple3(mesh.Force.coeff(i * 2),
                                mesh.Force.coeff(i * 2 + 1), 0);
    }
    pimpl->Grid->GetPointData()->AddArray(Force);

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
    pimpl->Grid->GetCellData()->AddArray(Normal);

    vtkNew<vtkFloatArray> Strain;
    Strain->SetNumberOfComponents(3);
    Strain->SetName("Strain");
    for (auto&& node : mesh.Nodes) {
        Strain->InsertNextTuple(node->Strain.data());
    }
    pimpl->Grid->GetPointData()->AddArray(Strain);

    vtkNew<vtkFloatArray> Stress;
    Stress->SetNumberOfComponents(3);
    Stress->SetName("Stress");
    for (auto&& node : mesh.Nodes) {
        Stress->InsertNextTuple(node->Stress.data());
    }
    pimpl->Grid->GetPointData()->AddArray(Stress);

    vtkNew<vtkFloatArray> Material;
    Material->SetNumberOfComponents(1);
    Material->SetName("Material");
    for (auto&& element : mesh.Elements) {
        Material->InsertNextValue(element->material->getIndex());
    }
    pimpl->Grid->GetCellData()->AddArray(Material);

    vtkNew<vtkFloatArray> AreaArray;
    AreaArray->SetNumberOfComponents(1);
    AreaArray->SetName("Area");
    for (auto&& element : mesh.Elements) {
        AreaArray->InsertNextValue(element->getArea());
    }
    pimpl->Grid->GetCellData()->AddArray(AreaArray);

    vtkNew<vtkFloatArray> strainEnergy;
    strainEnergy->SetNumberOfComponents(1);
    strainEnergy->SetName("StrainEnergy");
    for (auto&& element : mesh.Elements) {
        strainEnergy->InsertNextValue(element->getStrainEnergy());
    }
    pimpl->Grid->GetCellData()->AddArray(strainEnergy);

    vtkNew<vtkIntArray> Boundary;
    Boundary->SetNumberOfComponents(1);
    Boundary->SetName("Boundary");
    for (auto&& node : mesh.Nodes) {
        Boundary->InsertNextValue(node->isBoundary);
    }
    pimpl->Grid->GetPointData()->AddArray(Boundary);
}

void vtkManager::write(std::string fileName, bool isBinary) const {
    pimpl->writer->SetFileName(fileName.c_str());
    pimpl->writer->SetInputData(pimpl->Grid);
    if (isBinary) {
        pimpl->writer->SetDataModeToBinary();
    } else {
        pimpl->writer->SetDataModeToAscii();
    }
    pimpl->writer->Write();
}