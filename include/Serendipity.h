#pragma once
#include <memory>
#include <vector>

#include "Element.h"

class Serendipity : public Element {
   public:
    static const size_t nodeNum = 8;
    bool planeStress;

    explicit Serendipity(const size_t index_)
        : Element(index_, "Serendipity"), planeStress(true) {
        nodes = std::vector<std::shared_ptr<Node>>(nodeNum, nullptr);
    };
    Serendipity(size_t index_, const std::vector<std::shared_ptr<Node>>& nodes_,
                bool planeStress_ = true);
    ~Serendipity() override {};

    double getArea() const override;

    std::tuple<Eigen::VectorXd, Eigen::VectorXd> getShapeFuncLocalDerivative(
        double ksi, double eta) const override;

    std::tuple<Eigen::VectorXd, Eigen::VectorXd> getShapeFuncDerivative(
        double ksi, double eta) const override;

    Eigen::MatrixXd getJacobian(double ksi, double eta) const override;

    Eigen::MatrixXd getStrainMatrix(double ksi, double eta) const override;

    Eigen::MatrixXd getElasticMatrix() const override;

    Eigen::MatrixXd getStiffnessMatrix() const override;

    std::vector<Eigen::VectorXd> getGaussPointsStrain() const override;

    std::tuple<std::vector<Eigen::VectorXd>, std::vector<Eigen::VectorXd>>
    getGaussPointsStrainStress() const override;

    double getStrainEnergy() const override;

    void calculateNodeStrainStress() const override;
    void calculateNodeStrainStressViaGaussPoint() const override;
};