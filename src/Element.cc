#include "Element.h"
#include "Material.h"

Element::Element(const size_t index_, const std::string elementName_)
    : index(index_),
      elementName(elementName_),
      thickness(1.0),
      material(nullptr){};

int Element::setMaterial(Material* material) {
    if (material == nullptr) {
        return 1;
    }
    this->material = material;

    return 0;
}