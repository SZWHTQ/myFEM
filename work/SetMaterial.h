#include <vector>

#include "Element.h"
#include "Material.h"
#include "Mesh.h"

std::vector<double> get_element_center(Element* element);

void set_material(Mesh* mesh, std::vector<Material*> materials, double a,
                  double b);