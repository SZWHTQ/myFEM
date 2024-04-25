#pragma once
#include <memory>
#include <vector>

#include "Boundary.h"
#include "Element.h"

typedef std::pair<size_t, size_t> Key;

struct KeyHash {
    std::size_t operator()(const Key& k) const {
        return std::hash<size_t>()(k.first) ^ std::hash<size_t>()(k.second);
    }
};

struct KeyEqual {
    bool operator()(const Key& lhs, const Key& rhs) const {
        return lhs.first == rhs.first && lhs.second == rhs.second;
    }
};

class Mesh {
   public:
    size_t nodeNum;
    std::vector<std::shared_ptr<Node>> Nodes;
    std::vector<Element*> Elements;

    bool planeStress;
    enum MeshType { null, serendipity };

    Eigen::SparseVector<double> Force;

    Mesh() : nodeNum(0), planeStress(true) {
        Nodes = std::vector<std::shared_ptr<Node>>();
        Elements = std::vector<Element*>();
    }
    Mesh(MeshType type, std::vector<double> nodeCoord,
         std::vector<size_t> elementNodeTags,
         std::vector<size_t> boundaryNodeTags, bool planeStress = true);
    Mesh(std::vector<double> nodeCoord,
         std::vector<std::pair<MeshType, std::vector<size_t>>>
             elementTypeAndNodeTags,
         std::vector<size_t> boundaryNodeTags, bool planeStress = true);
    ~Mesh();

    Eigen::MatrixXd assembleStiffnessMatrix();
    Eigen::SparseMatrix<double> sparseAssembleStiffnessMatrix();
    Eigen::MatrixXd parallelAssembleStiffnessMatrix();
    Eigen::SparseMatrix<double> parallelSparseAssembleStiffnessMatrix();
    std::unordered_map<Key, double, KeyHash, KeyEqual> getStiffnessMatrixMap();
    static std::vector<double> const equivalentForce(Load* load);

    int Solve(std::list<Load>& loads, std::list<Boundary>& boundaries,
              bool verbose = false);

    int cuSolver(std::list<Load>& loads, std::list<Boundary>& boundaries,
                 bool verbose = false);
};
