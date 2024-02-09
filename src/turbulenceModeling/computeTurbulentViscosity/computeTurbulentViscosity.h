//
// Created by ruben on 8/02/24.
//

#ifndef SOLVERNS_COMPUTETURBULENTVISCOSITY_H
#define SOLVERNS_COMPUTETURBULENTVISCOSITY_H

#include "../../fields/scalarField/ScalarField.h"
#include "../../fields/vectorField/VectorField.h"
#include "../LES/modelLES.h"


ScalarField computeTurbulentViscosity(const PolyMesh& theMesh, const VectorField& u, const VectorBoundaryConditions& uBCs, modelLES turbulenceModel);

#endif //SOLVERNS_COMPUTETURBULENTVISCOSITY_H