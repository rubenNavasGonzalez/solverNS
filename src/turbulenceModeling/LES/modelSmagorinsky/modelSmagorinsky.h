//
// Created by ruben on 8/02/24.
//

#ifndef SOLVERNS_MODELSMAGORINSKY_H
#define SOLVERNS_MODELSMAGORINSKY_H

#include "../../computeTurbulentViscosity/computeTurbulentViscosity.h"


ScalarField modelSmagorinsky(const PolyMesh& theMesh, const VectorField& u, const VectorBoundaryConditions& uBCs);

#endif //SOLVERNS_MODELSMAGORINSKY_H
