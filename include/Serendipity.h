#pragma once
#include <memory>
#include <vector>

#include "Element.h"

class Serendipity : public Element {
   public:
    static size_t nodeNum;
    bool planeStress;

    explicit Serendipity(const size_t index_)
        : Element(index_, "Serendipity"), planeStress(true) {
        nodes = std::vector<std::shared_ptr<Node>>(nodeNum, nullptr);
    };
    Serendipity(size_t index_, const std::vector<std::shared_ptr<Node>>& nodes_,
                bool planeStress_ = true);

    [[nodiscard]] double getArea() const override;

    [[nodiscard]] const std::tuple<Eigen::VectorXd, Eigen::VectorXd> getShapeFuncLocalDerivative(
        double ksi, double eta) const override;

    [[nodiscard]] const std::tuple<Eigen::VectorXd, Eigen::VectorXd> getShapeFuncDerivative(
        double ksi, double eta) const override;

    [[nodiscard]] const Eigen::MatrixXd getJacobian(double ksi, double eta) const override;

    [[nodiscard]] const Eigen::MatrixXd getStrainMatrix(double ksi, double eta) const override;

    [[nodiscard]] const Eigen::MatrixXd getElasticMatrix() const override;

    [[nodiscard]] const Eigen::MatrixXd getStiffnessMatrix() const override;

    [[nodiscard]] const std::vector<Eigen::VectorXd> getGaussPointsStrain() const override;

    [[nodiscard]] const std::vector<Eigen::VectorXd> getGaussPointsStress() const override;

    [[nodiscard]] double getStrainEnergy() const override;

    void calculateStrainStress() const override;
    void calculateStrainStressGaussPoint() const override;

    ~Serendipity() override{};
};