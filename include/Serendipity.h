#pragma once
#include <memory>
#include <vector>

#include "Element.h"

class Serendipity : public Element {
   public:
    static size_t nodeNum;
    bool planeStress;

    Serendipity(const size_t index_)
        : Element(index_, "Serendipity"), planeStress(true) {
        nodes = std::vector<std::shared_ptr<Node>>(nodeNum, nullptr);
    };
    Serendipity(const size_t index_, const std::vector<std::shared_ptr<Node>>& nodes_,
                const bool planeStress_ = true);

    double getArea() const override;

    const std::tuple<Eigen::VectorXd, Eigen::VectorXd> getShapeFuncLocalDerivative(
        double ksi, double eta) const override;

    const std::tuple<Eigen::VectorXd, Eigen::VectorXd> getShapeFuncDerivative(
        double ksi, double eta) const override;

    const Eigen::MatrixXd getJacobian(double ksi, double eta) const override;

    const Eigen::MatrixXd getStrainMatrix(double ksi, double eta) const override;

    const Eigen::MatrixXd getElasticMatrix() const override;

    const Eigen::MatrixXd getStiffnessMatrix() const override;

    const std::vector<Eigen::VectorXd> getGaussPointsStrain() const override;

    const std::vector<Eigen::VectorXd> getGaussPointsStress() const override;

    double getStrainEnergy() const override;

    int calculateStrainStress() const override;
    int calculateStrainStressGaussPoint() const override;

    ~Serendipity(){};
};