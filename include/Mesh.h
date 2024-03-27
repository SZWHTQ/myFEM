#pragma once
#include <vector>

#include "Element.h"
#include "Boundary.h"

class Mesh {
   public:
    std::vector<Node> Nodes;
    std::vector<Element*> Elements;

    bool planeStress;
    enum MeshType { serendipity } meshType;

    Mesh() {
        Nodes = std::vector<Node>();
        Elements = std::vector<Element*>();
    }
    Mesh(MeshType type, std::vector<double> nodeCoord,
         std::vector<size_t> elementNodeTags, bool planeStress = true);

    Eigen::MatrixXd assembleStiffnessMatrix();
    static const std::vector<double> equivalentForce(Load* load);

    int Solve(std::list<Load>& loads, std::list<Boundary>& boundaries);
};