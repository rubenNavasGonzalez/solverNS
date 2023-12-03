//
// Created by ruben on 26/11/23.
//

#ifndef FLOWBETWEENFLATPLATES_CONVECTIVEORTHOGONAL_H
#define FLOWBETWEENFLATPLATES_CONVECTIVEORTHOGONAL_H

#include "../../fields/scalarField/ScalarField.h"
#include "../../fields/vectorField/VectorField.h"
#include "../../mesh/polyMesh/PolyMesh.h"
#include "../../boundaryConditions/vectorBoundaryConditions/VectorBoundaryConditions.h"


namespace fvc {

    VectorField convectiveOrthogonal(const ScalarField& mDot, const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs, const std::string& scheme);
}

#endif //FLOWBETWEENFLATPLATES_CONVECTIVEORTHOGONAL_H
