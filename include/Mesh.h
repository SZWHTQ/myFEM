#pragma once
#include <memory>
#include <vector>

#include "Element.h"
#include "Boundary.h"

class Mesh {
   public:
    std::vector<std::shared_ptr<Node>> Nodes;
    std::vector<Element*> Elements;

    bool planeStress;
    enum MeshType { serendipity } meshType;

    Eigen::VectorXd Force;

    Mesh() {
        Nodes = std::vector<std::shared_ptr<Node>>();
        Elements = std::vector<Element*>();
    }
    Mesh(MeshType type, std::vector<double> nodeCoord,
         std::vector<size_t> elementNodeTags, bool planeStress = true);

    Eigen::MatrixXd assembleStiffnessMatrix();
    static const std::vector<double> equivalentForce(Load* load);

    int Solve(std::list<Load>& loads, std::list<Boundary>& boundaries);
};