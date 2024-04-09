#pragma once
#include <Eigen/Dense>
#include <memory>

#include "Element.h"
#include "Mesh.h"

class Line {
   public:
    std::shared_ptr<Node> Start;
    std::shared_ptr<Node> End;
    std::shared_ptr<Node> Middle;
    Eigen::Vector3d Normal;
};

void getStrainEnergy(Mesh* mesh);