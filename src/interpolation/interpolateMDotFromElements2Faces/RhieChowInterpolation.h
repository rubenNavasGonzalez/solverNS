//
// Created by ruben on 11/12/23.
//

#ifndef FLOWBETWEENFLATPLATES_RHIECHOWINTERPOLATION_H
#define FLOWBETWEENFLATPLATES_RHIECHOWINTERPOLATION_H

#include "../../fields/scalarField/ScalarField.h"
#include "../../fields/vectorField/VectorField.h"
#include "../../boundaryConditions/scalarBoundaryConditions/ScalarBoundaryConditions.h"
#include "../../boundaryConditions/vectorBoundaryConditions/VectorBoundaryConditions.h"


ScalarField RhieChowInterpolation(const VectorField& u, const ScalarField& p, const PolyMesh& theMesh, double DeltaT, const VectorBoundaryConditions& uBCs, const ScalarBoundaryConditions& pBCs);

#endif //FLOWBETWEENFLATPLATES_RHIECHOWINTERPOLATION_H
