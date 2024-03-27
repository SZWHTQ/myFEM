#pragma once
#include <memory>

#include "Element.h"

class Serendipity : public Element {
   public:
    static const size_t nodeNum = 8;
    bool planeStress;

    Serendipity(const size_t index_)
        : Element(index_, "Serendipity"), planeStress(true) {
        nodes = std::vector<std::shared_ptr<Node>>(nodeNum, nullptr);
    };
    Serendipity(const size_t index_, const std::vector<std::shared_ptr<Node>>& nodes_,
                const bool planeStress_ = true);

    const std::tuple<Eigen::VectorXd, Eigen::VectorXd> shapeFuncLocalDerivative(
        double ksi, double eta) override;

    const std::tuple<Eigen::VectorXd, Eigen::VectorXd> shapeFuncDerivative(
        double ksi, double eta) override;

    const Eigen::MatrixXd Jacobian(double ksi, double eta) override;

    const Eigen::MatrixXd strainMatrix(double ksi, double eta) override;

    const Eigen::MatrixXd elasticMatrix(bool planeStress = true) override;

    const Eigen::MatrixXd stiffnessMatrix() override;

    ~Serendipity(){};
};