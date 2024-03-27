#pragma once
#include "Boundary.h"
#include "Mesh.h"

std::list<Load> apply_load(Mesh* mesh, double L, double B);

std::list<Boundary> apply_boundary(Mesh* mesh, double X, double Y);