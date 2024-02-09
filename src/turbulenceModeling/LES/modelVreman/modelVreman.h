//
// Created by ruben on 8/02/24.
//

#ifndef SOLVERNS_MODELVREMAN_H
#define SOLVERNS_MODELVREMAN_H

#include "../../computeTurbulentViscosity/computeTurbulentViscosity.h"


ScalarField modelVreman(const PolyMesh& theMesh, const VectorField& u, const VectorBoundaryConditions& uBCs);

#endif //SOLVERNS_MODELVREMAN_H
