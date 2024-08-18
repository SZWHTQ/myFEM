#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>

#include "Element.h"
#include "GaussIntegral.h"
#include "Material.h"
#include "Serendipity.h"

const size_t Serendipity::nodeNum;
using GaussIntegral::GaussData;

static constexpr std::array<int, 9> ksiId = {0, 0, 2, 2, 0, 1, 2, 1, 1};
static constexpr std::array<int, 9> etaId = {0, 2, 2, 0, 1, 2, 1, 0, 1};

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
    constexpr int gaussPointNum = 3;
    for (size_t i = 0; i < ksiId.size(); ++i) {
        auto&& J = getJacobian(GaussData<gaussPointNum>::abscissas[ksiId[i]],
                               GaussData<gaussPointNum>::abscissas[etaId[i]]);
        double detJ = J.determinant();
        area += detJ * GaussData<gaussPointNum>::weights[ksiId[i]] *
                GaussData<gaussPointNum>::weights[etaId[i]];
    }
    return area;
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd>
Serendipity::getShapeFuncLocalDerivative(const double ksi, const double eta) const {
    Eigen::Vector<double, nodeNum> shapeFunction_ksi, shapeFunction_eta;
    constexpr std::array<int, 4> k = {-1, 1, 1, -1}, e = {-1, -1, 1, 1};

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
Serendipity::getShapeFuncDerivative(const double ksi, const double eta) const {
    const auto J = getJacobian(ksi, eta);
    const double detJ = J.determinant();
    const auto [N_ksi, N_eta] = getShapeFuncLocalDerivative(ksi, eta);
    Eigen::VectorXd shapeFunction_x =
        (J(1, 1) * N_ksi - J(1, 0) * N_eta) / detJ;
    Eigen::VectorXd shapeFunction_y =
        (-J(0, 1) * N_ksi + J(0, 0) * N_eta) / detJ;

    return {shapeFunction_x, shapeFunction_y};
}

Eigen::MatrixXd Serendipity::getJacobian(const double ksi, const double eta) const {
    Eigen::Matrix<double, 2, 2> J;
    J.setZero();
    const auto [N_ksi, N_eta] = getShapeFuncLocalDerivative(ksi, eta);
    for (size_t i = 0; i < nodeNum; ++i) {
        J(0, 0) += N_ksi(i) * nodes.at(i)->x;
        J(1, 0) += N_ksi(i) * nodes.at(i)->y;
        J(0, 1) += N_eta(i) * nodes.at(i)->x;
        J(1, 1) += N_eta(i) * nodes.at(i)->y;
    }
    return J;
}

Eigen::MatrixXd Serendipity::getStrainMatrix(const double ksi, const double eta) const {
    const auto [N_x, N_y] = getShapeFuncDerivative(ksi, eta);
    Eigen::Matrix<double, 3, nodeNum * 2> strainMatrix;
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
    const double E = material->getElasticModulus();
    const double nu = material->getPoissonRatio();

    double A1, A2, A3;
    Eigen::Matrix<double, 3, 3> elasticMatrix;
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
    Eigen::Matrix<double, nodeNum * 2, nodeNum * 2> stiffnessMatrix;
    stiffnessMatrix.setZero();
    const auto D = getElasticMatrix();

    constexpr int gaussPointNum = 3;
    for (int i = 0; i < gaussPointNum; ++i) {
        for (int j = 0; j < gaussPointNum; ++j) {
            const auto B =
                getStrainMatrix(GaussData<gaussPointNum>::abscissas[i],
                                GaussData<gaussPointNum>::abscissas[j]);
            const double detJ =
                getJacobian(GaussData<gaussPointNum>::abscissas[i],
                            GaussData<gaussPointNum>::abscissas[j])
                    .determinant();
            stiffnessMatrix += thickness * B.transpose() * D * B * detJ *
                               GaussData<gaussPointNum>::weights[i] *
                               GaussData<gaussPointNum>::weights[j];
        }
    }

    return stiffnessMatrix;
}

void Serendipity::calculateNodeStrainStress() const {
    const auto D = getElasticMatrix();
    constexpr std::array<int, nodeNum> Ksi = {-1, 1, 1, -1, 0, 1, 0, -1};
    constexpr std::array<int, nodeNum> Eta = {-1, -1, 1, 1, -1, 0, 1, 0};
    Eigen::VectorXd displacementArray(nodeNum * 2);
    for (size_t i = 0; i < nodeNum; ++i) {
        displacementArray(2 * i) = nodes[i]->Displacement(0);
        displacementArray(2 * i + 1) = nodes[i]->Displacement(1);
    }
    for (size_t n = 0; n < nodeNum; ++n) {
        const auto B = getStrainMatrix(Ksi[n], Eta[n]);
        nodes[n]->Strain = B * displacementArray;
        nodes[n]->Stress = D * nodes[n]->Strain;
    }
}

std::vector<Eigen::VectorXd> Serendipity::getGaussPointsStrain() const {
    Eigen::Vector<double, nodeNum * 2> displacementArray;
    for (size_t i = 0; i < nodeNum; ++i) {
        displacementArray(2 * i) = nodes[i]->Displacement(0);
        displacementArray(2 * i + 1) = nodes[i]->Displacement(1);
    }

    // Calculate strain at gauss point
    constexpr int gaussPointNum = 3;
    std::vector<Eigen::VectorXd> gaussStrain(gaussPointNum * gaussPointNum,
                                             Eigen::VectorXd::Zero(3));

    for (size_t i = 0; i < ksiId.size(); ++i) {
        gaussStrain[i] =
            getStrainMatrix(GaussData<gaussPointNum>::abscissas[ksiId[i]],
                            GaussData<gaussPointNum>::abscissas[etaId[i]]) *
            displacementArray;
    }

    return gaussStrain;
}

std::tuple<std::vector<Eigen::VectorXd>, std::vector<Eigen::VectorXd>>
Serendipity::getGaussPointsStrainStress() const {
    const auto D = getElasticMatrix();
    const auto strainGp = getGaussPointsStrain();

    std::vector<Eigen::VectorXd> stressGp(strainGp.size());
    for (size_t i = 0; i < strainGp.size(); ++i) {
        stressGp[i] = D * strainGp[i];
    }

    return {strainGp, stressGp};
}

double Serendipity::getStrainEnergy() const {
    // auto& strainGp = getGaussPointsStrain();
    const auto [strainGp, stressGp] = getGaussPointsStrainStress();
    constexpr int gaussPointNum = 3;

    double strainEnergy = 0;
    for (size_t i = 0; i < ksiId.size(); ++i) {
        const auto J =
            getJacobian(GaussData<gaussPointNum>::abscissas[ksiId[i]],
                        GaussData<gaussPointNum>::abscissas[etaId[i]]);
        const double detJ = J.determinant();
        strainEnergy += 0.5 * stressGp[i].dot(strainGp[i]) * detJ *
                        GaussData<gaussPointNum>::weights[ksiId[i]] *
                        GaussData<gaussPointNum>::weights[etaId[i]];
    }

    return strainEnergy;
}

void Serendipity::calculateNodeStrainStressViaGaussPoint() const {
    Eigen::Vector<double, nodeNum * 2> displacementArray;
    for (size_t i = 0; i < nodeNum; ++i) {
        displacementArray(2 * i) = nodes[i]->Displacement(0);
        displacementArray(2 * i + 1) = nodes[i]->Displacement(1);
    }

    // Calculate strain at gauss point
    constexpr int gaussPointNum = 3;
    const auto strainGp = getGaussPointsStrain();

    // Calculate strain and stress at node
    const auto D = getElasticMatrix();
    Eigen::MatrixXd interpolationMatrix(nodeNum, nodeNum);
    interpolationMatrix.setZero();
    constexpr std::array<int, nodeNum> R = {-1, 1, 1, -1, 0, 1, 0, -1};
    constexpr std::array<int, nodeNum> S = {-1, -1, 1, 1, -1, 0, 1, 0};
    for (size_t n = 0; n < nodeNum; ++n) {
        const double r = R[n] / (-GaussData<gaussPointNum>::abscissas[0]),
                     s = S[n] / (-GaussData<gaussPointNum>::abscissas[0]);

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

    auto solver = interpolationMatrix.colPivHouseholderQr();
    Eigen::Vector<double, nodeNum> strain;
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
