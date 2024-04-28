#include "Element.h"

#include <utility>
#include "Material.h"

Element::Element(const size_t index_, std::string elementName_)
    : index(index_),
      elementName(std::move(elementName_)),
      material(nullptr),
      thickness(1.0){}

size_t Element::setMaterial(Material* material_) {
    if (material_ == nullptr) {
        return 1;
    }
    this->material = material_;

    return 0;
}