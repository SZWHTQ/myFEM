#pragma once
class Material {
   public:
    double youngModulus, poissonRatio;
    Material(double E_, double nu_) : youngModulus(E_), poissonRatio(nu_){};
    ~Material(){};
};

class Elastic : public Material {
   public:
    Elastic(double E_, double nu_) : Material(E_, nu_){};
    ~Elastic(){};
};