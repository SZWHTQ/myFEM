#pragma once
#include <vtkPoints.h>
#include <vtkQuadraticQuad.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>

#include <string>

#include "Mesh.h"

class vtkManager {
   public:
    vtkNew<vtkUnstructuredGrid> Grid;
    vtkNew<vtkUnstructuredGridWriter> writer;

    vtkManager(){};
    vtkManager(Mesh& mesh);
    ~vtkManager(){};

    void setData(std::string name, Mesh& mesh);
    void write(std::string fileName);
};