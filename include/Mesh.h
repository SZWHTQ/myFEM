#pragma once
#include <Eigen/Eigen>
#include <memory>
#include <vector>

class Load;
class Boundary;
class Node;
class Element;

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
         const std::vector<size_t>& boundaryNodeTags, bool planeStress = true);
    Mesh(std::vector<double> nodeCoord,
         std::vector<std::pair<MeshType, std::vector<size_t>>>
             elementsTypeAndNodeTags,
         const std::vector<size_t>& boundaryNodeTags, bool planeStress = true);
    ~Mesh();

    [[maybe_unused]] Eigen::MatrixXd const assembleStiffnessMatrix();
    Eigen::SparseMatrix<double> assembleSparseStiffnessMatrix();

    [[maybe_unused]] Eigen::MatrixXd parallelAssembleStiffnessMatrix();
    [[maybe_unused]] Eigen::SparseMatrix<double> parallelAssembleSparseStiffnessMatrix();
    static std::vector<double> const equivalentForce(Load* load);

    int Solve(std::list<Load>& loads, std::list<Boundary>& boundaries,
              bool verbose = false);
};

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& pair) const {
        auto hash1 = std::hash<T1>{}(pair.first);
        auto hash2 = std::hash<T2>{}(pair.second);
        return hash1 ^ hash2;
    }
};