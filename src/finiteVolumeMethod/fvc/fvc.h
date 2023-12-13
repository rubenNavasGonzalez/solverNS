//
// Created by ruben on 7/12/23.
//

#ifndef FLOWBETWEENFLATPLATES_FVC_H
#define FLOWBETWEENFLATPLATES_FVC_H

#include "../../fields/vectorField/VectorField.h"
#include "../../fields/scalarField/ScalarField.h"
#include "../../boundaryConditions/scalarBoundaryConditions/ScalarBoundaryConditions.h"


namespace fvc {

    VectorField convectiveOrthogonal(const ScalarField& mDot, const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs, const std::string& scheme);

    ScalarField divergence(const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs);

    VectorField gradient(const ScalarField& Phi, const PolyMesh& theMesh, const ScalarBoundaryConditions& PhiBCs);

    VectorField laplacianOrthogonal(const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs);
}

#endif //FLOWBETWEENFLATPLATES_FVC_H
