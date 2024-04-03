#include <iostream>
#include <memory>
#include <vector>

#include "Element.h"
#include "Material.h"
#include "Serendipity.h"

int main() {
    std::vector<std::shared_ptr<Node>> nodes;
    nodes.push_back(std::make_shared<Node>(0, 0, -1));
    nodes.push_back(std::make_shared<Node>(1, 4, 0));
    nodes.push_back(std::make_shared<Node>(2, 2, 3));
    nodes.push_back(std::make_shared<Node>(3, -2.5, 2.5));
    nodes.push_back(std::make_shared<Node>(4, 2, -0.5));
    nodes.push_back(std::make_shared<Node>(5, 3, 1.5));
    nodes.push_back(std::make_shared<Node>(6, -0.25, 2.75));
    nodes.push_back(std::make_shared<Node>(7, -1.25, 0.75));

    Elastic material(1, 0.3);

    std::shared_ptr<Element> counterClockWise =
        std::make_shared<Serendipity>(0, nodes);
    counterClockWise->setMaterial(&material);
    auto& K1 = counterClockWise->getStiffnessMatrix();
    std::cout << "Counter clockwise:" << std::endl;
    std::cout << K1 << std::endl;

    nodes = {nodes[3], nodes[2], nodes[1], nodes[0],
             nodes[7], nodes[6], nodes[5], nodes[4]};

    std::shared_ptr<Element> clockWise =
        std::make_shared<Serendipity>(1, nodes);
    clockWise->setMaterial(&material);
    auto& K2 = clockWise->getStiffnessMatrix();
    std::cout << "Clockwise:" << std::endl;
    std::cout << K2 << std::endl;

    std::cout << "Difference:" << std::endl;
    std::cout << K1 - K2 << std::endl;
    std::cout << "Norm: " << (K1 - K2).norm() << std::endl;
    std::cout << "Determinant: " << K1.determinant() << " " << K2.determinant()
              << std::endl;
}