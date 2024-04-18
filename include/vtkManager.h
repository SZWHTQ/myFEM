#pragma once

#include <string>

#include "Mesh.h"

// Declaration of an implementation class
class vtkManagerImpl;

class vtkManager {
   public:
    vtkManager();
    vtkManager(Mesh& mesh);
    ~vtkManager();

    void setMeshData(Mesh& mesh) const;
    void write(std::string fileName, bool isBinary = false) const;

   private:
    std::unique_ptr<vtkManagerImpl> pimpl;
};
