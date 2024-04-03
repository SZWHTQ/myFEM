#pragma once
#include <cstddef>
class Material {
   private:
    size_t index;
    double elasticModulus, poissonRatio;

   public:
    Material(size_t index_, double E_, double nu_)
        : index(index_), elasticModulus(E_), poissonRatio(nu_){};
    ~Material(){};
    size_t getIndex() { return index; };
    double getElasticModulus() { return elasticModulus; };
    double getPoissonRatio() { return poissonRatio; };
};

class Elastic : public Material {
   public:
    Elastic(size_t id, double E, double nu) : Material(id, E, nu){};
    ~Elastic(){};
};