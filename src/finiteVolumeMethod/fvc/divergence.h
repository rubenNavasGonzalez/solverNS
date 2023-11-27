//
// Created by ruben on 26/11/23.
//

#ifndef FLOWBETWEENFLATPLATES_DIVERGENCE_H
#define FLOWBETWEENFLATPLATES_DIVERGENCE_H

#include "../../fields/scalarField/ScalarField.h"
#include "../../fields/vectorField/VectorField.h"
#include "../../mesh/polyMesh/PolyMesh.h"
#include "../../boundaryConditions/vectorBoundaryConditions/VectorBoundaryConditions.h"


namespace fvc {

    ScalarField divergence(const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs);
}

#endif //FLOWBETWEENFLATPLATES_DIVERGENCE_H
