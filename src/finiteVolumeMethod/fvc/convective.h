//
// Created by ruben on 26/11/23.
//

#ifndef FLOWBETWEENFLATPLATES_CONVECTIVE_H
#define FLOWBETWEENFLATPLATES_CONVECTIVE_H

#include "../../fields/scalarField/ScalarField.h"
#include "../../fields/vectorField/VectorField.h"
#include "../../mesh/polyMesh/PolyMesh.h"
#include "../../boundaryConditions/vectorBoundaryConditions/VectorBoundaryConditions.h"


namespace fvc {

    VectorField convective(const ScalarField& mDot, const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs, const std::string& scheme);
}

#endif //FLOWBETWEENFLATPLATES_CONVECTIVE_H
