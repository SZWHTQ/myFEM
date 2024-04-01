#pragma once
#include <Eigen/Eigen>
#include <memory>
#include <string>
// #include <tuple>
#include <vector>

class Material;

class Node {
   private:
    size_t index;

   public:
    double x, y, z;

    Eigen::Vector3d Displacement{Eigen::Vector3d::Zero()};
    Eigen::Vector3d Strain{Eigen::Vector3d::Zero()};
    Eigen::Vector3d Stress{Eigen::Vector3d::Zero()};

    Node() : index(0), x(0), y(0), z(0){};
    Node(double x_, double y_, double z_) : x(x_), y(y_), z(z_){};
    Node(size_t index_, double x_, double y_, double z_)
        : index(index_), x(x_), y(y_), z(z_){};
    ~Node(){};
    int getIndex() { return index; };
    void setIndex(int i) { index = i; };
};

class Element {
   private:
    size_t index;
    std::string elementName;

   public:
    std::vector<std::shared_ptr<Node>> nodes;
    Material* material;
    double thickness;

    Element(){};
    explicit Element(const size_t index_, const std::string elementName_);
    virtual ~Element();

    // Derivative of shape function with respect to local coordinates
    virtual const std::tuple<Eigen::VectorXd, Eigen::VectorXd>
    shapeFuncLocalDerivative(double ksi, double eta) = 0;

    // Derivative of shape function with respect to global coordinates
    virtual const std::tuple<Eigen::VectorXd, Eigen::VectorXd>
    shapeFuncDerivative(double ksi, double eta) = 0;

    // Jacobian matrix
    virtual const Eigen::MatrixXd Jacobian(double ksi, double eta) = 0;

    // Strain matrix
    virtual const Eigen::MatrixXd strainMatrix(double ksi, double eta) = 0;

    // Elastic matrix
    virtual const Eigen::MatrixXd elasticMatrix(bool planeStress = true) = 0;

    // Stiffness matrix
    virtual const Eigen::MatrixXd stiffnessMatrix() = 0;

    virtual int calculateStrainStress() = 0;
    virtual int calculateStrainStressGaussPoint() = 0;

    std::string getElementName() { return elementName; };
    int getIndex() { return index; };
    int setMaterial(Material* material_);
};
