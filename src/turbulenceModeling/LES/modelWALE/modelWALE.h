//
// Created by ruben on 8/02/24.
//

#ifndef SOLVERNS_MODELWALE_H
#define SOLVERNS_MODELWALE_H

#include "../../computeTurbulentViscosity/computeTurbulentViscosity.h"


ScalarField modelWALE(const PolyMesh& theMesh, const VectorField& u, const VectorBoundaryConditions& uBCs);

#endif //SOLVERNS_MODELWALE_H
