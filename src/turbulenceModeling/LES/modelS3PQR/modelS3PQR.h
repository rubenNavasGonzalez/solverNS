//
// Created by ruben on 9/02/24.
//

#ifndef SOLVERNS_MODELS3PQR_H
#define SOLVERNS_MODELS3PQR_H

#include "../../computeTurbulentViscosity/computeTurbulentViscosity.h"


ScalarField modelS3PQR(const PolyMesh& theMesh, const VectorField& u, const VectorBoundaryConditions& uBCs, modelLES turbulenceModel);

#endif //SOLVERNS_MODELS3PQR_H