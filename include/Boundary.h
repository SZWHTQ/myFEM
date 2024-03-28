#pragma once
#include <memory>
#include <vector>

#include "Element.h"
class Load {
   public:
    std::vector<std::shared_ptr<Node>> nodes{3, nullptr};
    std::vector<double> normalForce{0, 0, 0};
    std::vector<double> shearForce{0, 0, 0};
    double thickness = 0;
};

class Boundary {
   public:
    std::shared_ptr<Node> node;
    std::vector<double> displacement{0, 0, 0};
    std::vector<bool> fixed{false, false, false}; 
};