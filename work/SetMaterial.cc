#include <vector>

#include "SetMaterial.h"

std::vector<double> get_element_center(Element* element) {
    std::vector<double> center(3, 0.0);
    for (auto& node : element->nodes) {
        center[0] += node->x;
        center[1] += node->y;
        center[2] += node->z;
    }
    center[0] /= element->nodes.size();
    center[1] /= element->nodes.size();
    center[2] /= element->nodes.size();
    return center;
}

static bool inTheEllipse(double x, double y, double a, double b) {
    return x * x / a / a + y * y / b / b < 1;
}

void set_material(Mesh* mesh, std::vector<Material*> materials, double a,
                  double b) {
    for (auto&& element : mesh->Elements) {
        auto&& center = get_element_center(element);
        if (inTheEllipse(center[0], center[1], a, b)) {
            element->material = materials[1];
//             std::cout << "Element " << element->getIndex() << " is in the inclusion"
//                       << std::endl;
        } else {
            element->material = materials[0];
//             std::cout << "Element " << element->getIndex() << " is in the matrix"
//                       << std::endl;
        }
    }
}

void set_material(Mesh* mesh, std::vector<Material*> materials,
                  std::vector<size_t> elementMaterialTags) {
    for (size_t i = 0; i < mesh->Elements.size(); ++i) {
        mesh->Elements[i]->material = materials[elementMaterialTags[i]];
    }
}