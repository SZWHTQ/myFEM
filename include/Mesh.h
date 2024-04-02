#pragma once
#include <memory>
#include <vector>

#include "Element.h"
#include "Boundary.h"

#define VERBOSE 1

class Mesh {
   public:
    std::vector<std::shared_ptr<Node>> Nodes;
    std::vector<Element*> Elements;

    bool planeStress;
    enum MeshType { serendipity } meshType;

    Eigen::SparseVector<double> Force;

    Mesh() {
        Nodes = std::vector<std::shared_ptr<Node>>();
        Elements = std::vector<Element*>();
    }
    Mesh(MeshType type, std::vector<double> nodeCoord,
         std::vector<size_t> elementNodeTags, bool planeStress = true);

    Eigen::MatrixXd assembleStiffnessMatrix();
    Eigen::SparseMatrix<double> sparseAssembleStiffnessMatrix();
    Eigen::MatrixXd parallelAssembleStiffnessMatrix();
    static const std::vector<double> equivalentForce(Load* load);

    int Solve(std::list<Load>& loads, std::list<Boundary>& boundaries);
};

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &pair) const {
        auto hash1 = std::hash<T1>{}(pair.first);
        auto hash2 = std::hash<T2>{}(pair.second);
        return hash1 ^ hash2;
    }
};