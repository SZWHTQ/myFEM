#pragma once
#include <vtkPoints.h>
#include <vtkQuadraticQuad.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <string>

#include "Mesh.h"

class vtkManager {
   public:
    vtkNew<vtkUnstructuredGrid> Grid;
    vtkNew<vtkXMLUnstructuredGridWriter> writer;

    vtkManager(){};
    vtkManager(Mesh& mesh);
    ~vtkManager(){};

    void setData(Mesh& mesh);
    void write(std::string fileName, bool isBinary = false);
};