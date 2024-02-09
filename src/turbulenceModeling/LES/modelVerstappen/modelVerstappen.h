//
// Created by ruben on 8/02/24.
//

#ifndef SOLVERNS_MODELVERSTAPPEN_H
#define SOLVERNS_MODELVERSTAPPEN_H

#include "../../computeTurbulentViscosity/computeTurbulentViscosity.h"


ScalarField modelVerstappen(const PolyMesh& theMesh, const VectorField& u, const VectorBoundaryConditions& uBCs);

#endif //SOLVERNS_MODELVERSTAPPEN_H
