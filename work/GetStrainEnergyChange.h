#pragma once
#include <Eigen/Dense>

#include "Element.h"
#include "Material.h"
#include "Mesh.h"

double getStrainEnergyChange(Mesh* mesh, Material* matrix, Material* inclusion,
                             bool isPlaneStress);

double getStrainEnergyChange(Mesh* mesh, Mesh* meshNoInclusion,
                             Material* matrix, Material* inclusion,
                             bool isPlaneStress);
