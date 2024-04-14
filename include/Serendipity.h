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

    const double getArea() override;

    const std::tuple<Eigen::VectorXd, Eigen::VectorXd> getShapeFuncLocalDerivative(
        double ksi, double eta) override;

    const std::tuple<Eigen::VectorXd, Eigen::VectorXd> getShapeFuncDerivative(
        double ksi, double eta) override;

    const Eigen::MatrixXd getJacobian(double ksi, double eta) override;

    const Eigen::MatrixXd getStrainMatrix(double ksi, double eta) override;

    const Eigen::MatrixXd getElasticMatrix() override;

    const Eigen::MatrixXd getStiffnessMatrix() override;

    const std::vector<Eigen::VectorXd> getGaussPointsStrain() override;

    const std::vector<Eigen::VectorXd> getGaussPointsStress() override;

    const double getStrainEnergy() override;

    int calculateStrainStress() override;
    int calculateStrainStressGaussPoint() override;

    ~Serendipity(){};
};