#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>

#include "Element.h"
#include "GaussIntegral.h"
#include "Material.h"
#include "Serendipity.h"

const size_t Serendipity::nodeNum;

Serendipity::Serendipity(const size_t index_,
                         const std::vector<std::shared_ptr<Node>>& nodes_,
                         const bool planeStress_)
    : Element(index_, "Serendipity"), planeStress(planeStress_) {
    if (nodes_.size() == 8) {
        nodes = nodes_;
    } else {
        std::cerr << "Serendipity element must have 8 nodes" << std::endl;
    }
}

// get the area of the element via gauss integral
double Serendipity::getArea() const {
    double area = 0;
    int gaussPointNum = 3;
    const auto& gaussData = GaussIntegral::getGaussData(gaussPointNum);
    const std::vector<int> ksiId = {0, 0, 2, 2, 0, 1, 2, 1, 1};
    const std::vector<int> etaId = {0, 2, 2, 0, 1, 2, 1, 0, 1};
    for (size_t i = 0; i < ksiId.size(); ++i) {
        auto&& J = getJacobian(gaussData.abscissas[ksiId[i]],
                               gaussData.abscissas[etaId[i]]);
        double detJ = J.determinant();
        area +=
            detJ * gaussData.weights[ksiId[i]] * gaussData.weights[etaId[i]];
    }
    return area;
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd>
Serendipity::getShapeFuncLocalDerivative(double ksi, double eta) const {
    Eigen::Vector<double, nodeNum> shapeFunction_ksi, shapeFunction_eta;
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

std::tuple<Eigen::VectorXd, Eigen::VectorXd>
Serendipity::getShapeFuncDerivative(double ksi, double eta) const {
    auto&& J = getJacobian(ksi, eta);
    double detJ = J.determinant();
    auto&& [N_ksi, N_eta] = getShapeFuncLocalDerivative(ksi, eta);
    Eigen::VectorXd shapeFunction_x =
        (J(1, 1) * N_ksi - J(1, 0) * N_eta) / detJ;
    Eigen::VectorXd shapeFunction_y =
        (-J(0, 1) * N_ksi + J(0, 0) * N_eta) / detJ;

    return {shapeFunction_x, shapeFunction_y};
}

Eigen::MatrixXd Serendipity::getJacobian(double ksi, double eta) const {
    Eigen::MatrixXd J(2, 2);
    J.setZero();
    auto&& [N_ksi, N_eta] = getShapeFuncLocalDerivative(ksi, eta);
    for (size_t i = 0; i < nodeNum; ++i) {
        J(0, 0) += N_ksi(i) * nodes.at(i)->x;
        J(1, 0) += N_ksi(i) * nodes.at(i)->y;
        J(0, 1) += N_eta(i) * nodes.at(i)->x;
        J(1, 1) += N_eta(i) * nodes.at(i)->y;
    }
    return J;
}

Eigen::MatrixXd Serendipity::getStrainMatrix(double ksi, double eta) const {
    auto&& [N_x, N_y] = getShapeFuncDerivative(ksi, eta);
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

Eigen::MatrixXd Serendipity::getElasticMatrix() const {
    double E = material->getElasticModulus();
    double nu = material->getPoissonRatio();

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

Eigen::MatrixXd Serendipity::getStiffnessMatrix() const {
    Eigen::MatrixXd stiffnessMatrix(nodeNum * 2, nodeNum * 2);
    stiffnessMatrix.setZero();
    auto&& D = getElasticMatrix();

    int gaussPointNum = 3;
    auto& gaussData = GaussIntegral::getGaussData(gaussPointNum);

    for (int i = 0; i < gaussPointNum; ++i) {
        for (int j = 0; j < gaussPointNum; ++j) {
            auto&& B =
                getStrainMatrix(gaussData.abscissas[i], gaussData.abscissas[j]);
            double detJ =
                getJacobian(gaussData.abscissas[i], gaussData.abscissas[j])
                    .determinant();
            stiffnessMatrix += thickness * B.transpose() * D * B * detJ *
                               gaussData.weights[i] * gaussData.weights[j];
        }
    }

    return stiffnessMatrix;
}

void Serendipity::calculateNodeStrainStress() const {
    auto&& D = getElasticMatrix();
    const std::vector<int> Ksi = {-1, 1, 1, -1, 0, 1, 0, -1};
    const std::vector<int> Eta = {-1, -1, 1, 1, -1, 0, 1, 0};
    Eigen::VectorXd displacementArray(nodeNum * 2);
    for (size_t i = 0; i < nodeNum; ++i) {
        displacementArray(2 * i) = nodes[i]->Displacement(0);
        displacementArray(2 * i + 1) = nodes[i]->Displacement(1);
    }
    for (size_t n = 0; n < nodeNum; ++n) {
        auto&& B = getStrainMatrix(Ksi[n], Eta[n]);
        nodes[n]->Strain = B * displacementArray;
        nodes[n]->Stress = D * nodes[n]->Strain;
    }
}

/*
int Serendipity::calculateStrainStress() {
    Eigen::VectorXd displacementArray(nodeNum * 2);
    for (int i = 0; i < nodeNum; ++i) {
        displacementArray(2 * i) = nodes[i]->Displacement(0);
        displacementArray(2 * i + 1) = nodes[i]->Displacement(1);
    }

    // Calculate strain at gauss point
    int gaussPointNum = 3;
    const auto& gaussData = GaussIntegral::getGaussData(gaussPointNum);

    std::vector<Eigen::Vector3d> strainAtGaussPoints(
        gaussPointNum * gaussPointNum, Eigen::Vector3d::Zero());
    strainAtGaussPoints[0] =
        strainMatrix(gaussData.abscissas[0], gaussData.abscissas[0]) *
        displacementArray;
    strainAtGaussPoints[1] =
        strainMatrix(gaussData.abscissas[1], gaussData.abscissas[0]) *
        displacementArray;
    strainAtGaussPoints[2] =
        strainMatrix(gaussData.abscissas[1], gaussData.abscissas[1]) *
        displacementArray;
    strainAtGaussPoints[3] =
        strainMatrix(gaussData.abscissas[0], gaussData.abscissas[1]) *
        displacementArray;

    // Calculate strain and stress at node
    auto& D = elasticMatrix(planeStress);
    const std::vector<int> R = {-1, 1, 1, -1, 0, 1, 0, -1};
    const std::vector<int> S = {-1, -1, 1, 1, -1, 0, 1, 0};
    for (int n = 0; n < nodeNum; ++n) {
        nodes[n]->Strain.setZero();
        double r = R[n] / (-gaussData.abscissas[0]),
               s = S[n] / (-gaussData.abscissas[0]);
        std::vector<double> N(8);
        N[0] = 0.25 * (1 - r) * (1 - s);
        N[1] = 0.25 * (1 + r) * (1 - s);
        N[2] = 0.25 * (1 + r) * (1 + s);
        N[3] = 0.25 * (1 - r) * (1 + s);
        nodes[n]->Strain += N[0] * strainAtGaussPoints[0];
        nodes[n]->Strain += N[1] * strainAtGaussPoints[1];
        nodes[n]->Strain += N[2] * strainAtGaussPoints[2];
        nodes[n]->Strain += N[3] * strainAtGaussPoints[3];

        nodes[n]->Stress = D * nodes[n]->Strain;
    }
    return 0;
}
 */

std::vector<Eigen::VectorXd> Serendipity::getGaussPointsStrain() const {
    Eigen::VectorXd displacementArray(nodeNum * 2);
    for (size_t i = 0; i < nodeNum; ++i) {
        displacementArray(2 * i) = nodes[i]->Displacement(0);
        displacementArray(2 * i + 1) = nodes[i]->Displacement(1);
    }

    // Calculate strain at gauss point
    int gaussPointNum = 3;
    const auto& gaussData = GaussIntegral::getGaussData(gaussPointNum);
    std::vector<Eigen::VectorXd> gaussStrain(gaussPointNum * gaussPointNum,
                                             Eigen::VectorXd::Zero(3));

    const std::vector<int> ksiId = {0, 0, 2, 2, 0, 1, 2, 1, 1};
    const std::vector<int> etaId = {0, 2, 2, 0, 1, 2, 1, 0, 1};
    for (size_t i = 0; i < ksiId.size(); ++i) {
        gaussStrain[i] = getStrainMatrix(gaussData.abscissas[ksiId[i]],
                                         gaussData.abscissas[etaId[i]]) *
                         displacementArray;
    }

    return gaussStrain;
}

std::tuple<std::vector<Eigen::VectorXd>, std::vector<Eigen::VectorXd>>
Serendipity::getGaussPointsStrainStress() const {
    auto&& D = getElasticMatrix();
    auto&& strainGp = getGaussPointsStrain();

    std::vector<Eigen::VectorXd> stressGp(strainGp.size(),
                                          Eigen::VectorXd::Zero(3));
    for (size_t i = 0; i < strainGp.size(); ++i) {
        stressGp[i] = D * strainGp[i];
    }

    return {strainGp, stressGp};
}

double Serendipity::getStrainEnergy() const {
    // auto& strainGp = getGaussPointsStrain();
    auto&& [strainGp, stressGp] = getGaussPointsStrainStress();
    int gaussPointNum = 3;
    const auto& gaussData = GaussIntegral::getGaussData(gaussPointNum);
    const std::vector<int> ksiId = {0, 0, 2, 2, 0, 1, 2, 1, 1};
    const std::vector<int> etaId = {0, 2, 2, 0, 1, 2, 1, 0, 1};
    double strainEnergy = 0;
    for (size_t i = 0; i < ksiId.size(); ++i) {
        auto&& J = getJacobian(gaussData.abscissas[ksiId[i]],
                               gaussData.abscissas[etaId[i]]);
        double detJ = J.determinant();
        strainEnergy += 0.5 * stressGp[i].dot(strainGp[i]) * detJ *
                        gaussData.weights[ksiId[i]] *
                        gaussData.weights[etaId[i]];
    }

    return strainEnergy;
}

void Serendipity::calculateNodeStrainStressViaGaussPoint() const {
    Eigen::VectorXd displacementArray(nodeNum * 2);
    for (size_t i = 0; i < nodeNum; ++i) {
        displacementArray(2 * i) = nodes[i]->Displacement(0);
        displacementArray(2 * i + 1) = nodes[i]->Displacement(1);
    }

    // Calculate strain at gauss point
    int gaussPointNum = 3;
    const auto& gaussData = GaussIntegral::getGaussData(gaussPointNum);
    const auto& strainGp = getGaussPointsStrain();

    // Calculate strain and stress at node
    auto&& D = getElasticMatrix();
    Eigen::MatrixXd interpolationMatrix(nodeNum, nodeNum);
    interpolationMatrix.setZero();
    const std::vector<int> R = {-1, 1, 1, -1, 0, 1, 0, -1};
    const std::vector<int> S = {-1, -1, 1, 1, -1, 0, 1, 0};
    for (size_t n = 0; n < nodeNum; ++n) {
        double r = R[n] / (-gaussData.abscissas[0]),
               s = S[n] / (-gaussData.abscissas[0]);

        interpolationMatrix(n, 4) = 0.5 * (1 - r * r) * (1 - s);
        interpolationMatrix(n, 5) = 0.5 * (1 + r) * (1 - s * s);
        interpolationMatrix(n, 6) = 0.5 * (1 - r * r) * (1 + s);
        interpolationMatrix(n, 7) = 0.5 * (1 - r) * (1 - s * s);
        interpolationMatrix(n, 0) =
            0.25 * (1 - r) * (1 - s) -
            0.5 * (interpolationMatrix(n, 4) + interpolationMatrix(n, 7));
        interpolationMatrix(n, 1) =
            0.25 * (1 + r) * (1 - s) -
            0.5 * (interpolationMatrix(n, 4) + interpolationMatrix(n, 5));
        interpolationMatrix(n, 2) =
            0.25 * (1 + r) * (1 + s) -
            0.5 * (interpolationMatrix(n, 5) + interpolationMatrix(n, 6));
        interpolationMatrix(n, 3) =
            0.25 * (1 - r) * (1 + s) -
            0.5 * (interpolationMatrix(n, 6) + interpolationMatrix(n, 7));
    }

    auto& solver = interpolationMatrix.colPivHouseholderQr();
    Eigen::VectorXd strain(nodeNum);
    for (int d = 0; d < 3; ++d) {
        for (size_t i = 0; i < nodeNum; ++i) {
            strain(i) = strainGp[i](d);
        }
        Eigen::VectorXd strainAtNode = solver.solve(strain);
        for (size_t i = 0; i < nodeNum; ++i) {
            nodes[i]->Strain(d) = strainAtNode(i);
        }
    }

    for (size_t i = 0; i < nodeNum; ++i) {
        nodes[i]->Stress = D * nodes[i]->Strain;
    }
}
