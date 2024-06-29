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

    const bool planeStress;
    enum MeshType { null, serendipity };

    Eigen::SparseVector<double> Force;

    Mesh() : nodeNum(0), planeStress(true) {
        Nodes = std::vector<std::shared_ptr<Node>>();
        Elements = std::vector<Element*>();
    }
    Mesh(const MeshType type, const std::vector<double> nodeCoord,
         const std::vector<size_t> elementNodeTags,
         const std::vector<size_t>& boundaryNodeTags,
         const bool planeStress = true);
    Mesh(const std::vector<double> nodeCoord,
         const std::vector<std::pair<MeshType, std::vector<size_t>>>
             elementsTypeAndNodeTags,
         const std::vector<size_t>& boundaryNodeTags,
         const bool planeStress = true);
    ~Mesh();

    [[maybe_unused]] Eigen::MatrixXd assembleStiffnessMatrix();
    Eigen::SparseMatrix<double> assembleSparseStiffnessMatrix();

    [[maybe_unused]] Eigen::MatrixXd parallelAssembleStiffnessMatrix();
    [[maybe_unused]] Eigen::SparseMatrix<double>
    parallelAssembleSparseStiffnessMatrix();
    static const std::vector<double> equivalentForce(const Load& load);

    int Solve(const std::list<Load>& loads,
              const std::list<Boundary>& boundaries,
              const bool verbose = false);
};

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& pair) const {
        auto hash1 = std::hash<T1>{}(pair.first);
        auto hash2 = std::hash<T2>{}(pair.second);
        return hash1 ^ hash2;
    }
};