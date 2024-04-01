#pragma once
class Material {
   private:
    double elasticModulus, poissonRatio;

   public:
    Material(double E_, double nu_) : elasticModulus(E_), poissonRatio(nu_){};
    ~Material(){};
    double getElasticModulus() { return elasticModulus; };
    double getPoissonRatio() { return poissonRatio; };
};

class Elastic : public Material {
   public:
    Elastic(double E_, double nu_) : Material(E_, nu_){};
    ~Elastic(){};
};