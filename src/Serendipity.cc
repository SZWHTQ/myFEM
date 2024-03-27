#include <iostream>
#include <memory>
#include <vector>

#include "Element.h"
#include "GaussIntegral.h"
#include "Material.h"
#include "Serendipity.h"

Serendipity::Serendipity(const size_t index_, const std::vector<std::shared_ptr<Node>>& nodes_,
                         const bool planeStress_)
    : Element(index_, "Serendipity"), planeStress(planeStress_) {
    if (nodes_.size() == 8) {
        nodes = nodes_;
    } else {
        std::cerr << "Serendipity element must have 8 nodes" << std::endl;
    }
}

const std::tuple<Eigen::VectorXd, Eigen::VectorXd>
Serendipity::shapeFuncLocalDerivative(double ksi, double eta) {
    Eigen::VectorXd shapeFunction_ksi(8), shapeFunction_eta(8);
    const std::vector<int> k{-1, 1, 1, -1}, e{-1, -1, 1, 1};

    for (size_t i = 0; i < 4; ++i) {
        shapeFunction_ksi(i) =
            k[i] * (1 + e[i] * eta) * (k[i] * ksi + e[i] * eta - 1) * 0.25 +
            k[i] * (1 + k[i] * ksi) * (1 + e[i] * eta) * 0.25;
        shapeFunction_eta(i) =
            e[i] * (1 + k[i] * ksi) * (k[i] * ksi + e[i] * eta - 1) * 0.25 +
            e[i] * (1 + k[i] * ksi) * (1 + e[i] * eta) * 0.25;
    }

    shapeFunction_ksi(4) = ksi * (eta - 1);
    shapeFunction_ksi(5) = 0.5 * (1 - eta * eta);
    shapeFunction_ksi(6) = -ksi * (eta + 1);
    shapeFunction_ksi(7) = 0.5 * (eta * eta - 1);

    shapeFunction_eta(4) = 0.5 * (ksi * ksi - 1);
    shapeFunction_eta(5) = -eta * (ksi + 1);
    shapeFunction_eta(6) = 0.5 * (1 - ksi * ksi);
    shapeFunction_eta(7) = eta * (ksi - 1);

    return {shapeFunction_ksi, shapeFunction_eta};
}
/*
const std::tuple<Eigen::VectorXd, Eigen::VectorXd>
Serendipity::shapeFuncDerivative(double ksi, double eta) {
    auto& J = Jacobian(ksi, eta);
    auto& invJ = J.inverse();
    auto& [N_ksi, N_eta] = shapeFuncLocalDerivative(ksi, eta);
    Eigen::VectorXd shapeFunction_x = invJ(0, 0) * N_ksi + invJ(1, 0) *
    N_eta; Eigen::VectorXd shapeFunction_y = invJ(0, 1) * N_ksi + invJ(1, 1)
    * N_eta;

    return {shapeFunction_x, shapeFunction_y};
}
 */
const std::tuple<Eigen::VectorXd, Eigen::VectorXd>
Serendipity::shapeFuncDerivative(double ksi, double eta) {
    auto& J = Jacobian(ksi, eta);
    double detJ = J.determinant();
    auto& [N_ksi, N_eta] = shapeFuncLocalDerivative(ksi, eta);
    Eigen::VectorXd shapeFunction_x =
        (J(1, 1) * N_ksi - J(1, 0) * N_eta) / detJ;
    Eigen::VectorXd shapeFunction_y =
        (-J(0, 1) * N_ksi + J(0, 0) * N_eta) / detJ;

    return {shapeFunction_x, shapeFunction_y};
}

const Eigen::MatrixXd Serendipity::Jacobian(double ksi, double eta) {
    Eigen::MatrixXd J(2, 2);
    J.setZero();
    auto& [N_ksi, N_eta] = shapeFuncLocalDerivative(ksi, eta);
    for (size_t i = 0; i < nodeNum; ++i) {
        J(0, 0) += N_ksi(i) * nodes.at(i)->x;
        J(1, 0) += N_ksi(i) * nodes.at(i)->y;
        J(0, 1) += N_eta(i) * nodes.at(i)->x;
        J(1, 1) += N_eta(i) * nodes.at(i)->y;
    }
    return J;
}

const Eigen::MatrixXd Serendipity::strainMatrix(double ksi, double eta) {
    auto& [N_x, N_y] = shapeFuncDerivative(ksi, eta);
    Eigen::MatrixXd strainMatrix(3, nodeNum * 2);
    strainMatrix.setZero();

    for (size_t i = 0; i < nodeNum; ++i) {
        strainMatrix(0, 2 * i) = N_x(i);
        strainMatrix(1, 2 * i + 1) = N_y(i);
        strainMatrix(2, 2 * i) = N_y(i);
        strainMatrix(2, 2 * i + 1) = N_x(i);
    }

    return strainMatrix;
}

const Eigen::MatrixXd Serendipity::elasticMatrix(bool planeStress) {
    double E = material->youngModulus;
    double nu = material->poissonRatio;

    double A1, A2, A3;
    Eigen::MatrixXd elasticMatrix(3, 3);
    if (planeStress) {
        A1 = nu;
        A2 = (1 - nu) * 0.5;
        A3 = E / (1 - nu * nu);
    } else {
        A1 = nu / (1 - nu);
        A2 = (1 - 2 * nu) / 2 / (1 - nu);
        A3 = E * (1 - nu) / (1 + nu) / (1 - 2 * nu);
    }
    elasticMatrix << A3, A3 * A1, 0, A3 * A1, A3, 0, 0, 0, A3 * A2;

    return elasticMatrix;
}

const Eigen::MatrixXd Serendipity::stiffnessMatrix() {
    Eigen::MatrixXd stiffnessMatrix(nodeNum * 2, nodeNum * 2);
    stiffnessMatrix.setZero();
    auto& D = elasticMatrix();

    int gaussPointNum = 3;
    auto& gauss = GaussIntegral::getGaussData(gaussPointNum);

    for (int i = 0; i < gaussPointNum; ++i) {
        for (int j = 0; j < gaussPointNum; ++j) {
            auto& B = strainMatrix(gauss.abscissas[i], gauss.abscissas[j]);
            double detJ =
                Jacobian(gauss.abscissas[i], gauss.abscissas[j]).determinant();
            stiffnessMatrix += B.transpose() * D * B * detJ * gauss.weights[i] *
                               gauss.weights[j];
        }
    }

    return stiffnessMatrix;
}