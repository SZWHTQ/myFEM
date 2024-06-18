#pragma once

#include <memory>
#include <string>

class Mesh;

class vtkManager {
   public:
    vtkManager() = delete;
    vtkManager(Mesh& mesh);
    ~vtkManager();

    void setMeshData(Mesh& mesh) const;
    void write(std::string fileName, bool isBinary = false) const;

   private:
    // Declaration of an implementation class
    class vtkManagerImpl;

    std::unique_ptr<vtkManagerImpl> pimpl;
};
