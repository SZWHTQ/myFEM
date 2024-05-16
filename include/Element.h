#pragma once
#include <Eigen/Eigen>
#include <cstddef>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

class Material;

class Node {
   private:
    size_t index;

   public:
    double x, y, z;
    bool isBoundary = false;

    Eigen::Vector3d Displacement{Eigen::Vector3d::Zero()};
    Eigen::Vector3d Strain{Eigen::Vector3d::Zero()};
    Eigen::Vector3d Stress{Eigen::Vector3d::Zero()};

    Node() : index(0), x(0), y(0), z(0) {};
    Node(double x_, double y_, double z_) : index(0), x(x_), y(y_), z(z_) {};
    Node(size_t index_, double x_, double y_, double z_)
        : index(index_), x(x_), y(y_), z(z_) {};
    ~Node() = default;
    [[nodiscard]] size_t getIndex() const { return index; };
    void setIndex(int i) { index = i; };
};

class Element {
   public:
    std::vector<std::shared_ptr<Node>> nodes;
    Material* material;
    double thickness;
    static const size_t nodeNum = 0;

    Element() : material(nullptr), thickness(1.0), index(0) {};
    explicit Element(size_t index_, std::string elementName_);
    virtual ~Element() = default;

    // Area of element
    virtual double getArea() const = 0;

    // Derivative of shape function with respect to local coordinates
    virtual std::tuple<Eigen::VectorXd, Eigen::VectorXd>
    getShapeFuncLocalDerivative(double ksi, double eta) const = 0;

    // Derivative of shape function with respect to global coordinates
    virtual std::tuple<Eigen::VectorXd, Eigen::VectorXd> getShapeFuncDerivative(
        double ksi, double eta) const = 0;

    // Jacobian matrix
    virtual Eigen::MatrixXd getJacobian(double ksi, double eta) const = 0;

    // Strain matrix
    virtual Eigen::MatrixXd getStrainMatrix(double ksi, double eta) const = 0;

    // Elastic matrix
    virtual Eigen::MatrixXd getElasticMatrix() const = 0;

    // Stiffness matrix
    virtual Eigen::MatrixXd getStiffnessMatrix() const = 0;

    // Gauss Strain
    virtual std::vector<Eigen::VectorXd> getGaussPointsStrain() const = 0;

    // Gauss Stress
    virtual std::tuple<std::vector<Eigen::VectorXd>,
                       std::vector<Eigen::VectorXd>>
    getGaussPointsStrainStress() const = 0;

    // Strain energy
    virtual double getStrainEnergy() const = 0;

    // Calculate strain and stress
    virtual void calculateNodeStrainStress() const = 0;
    virtual void calculateNodeStrainStressViaGaussPoint() const = 0;

    std::string getElementName() const { return elementName; };
    size_t getIndex() const { return index; };
    size_t setMaterial(Material* material_);

   private:
    size_t index;
    std::string elementName;
};
