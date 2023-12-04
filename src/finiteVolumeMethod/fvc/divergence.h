//
// Created by ruben on 4/12/23.
//

#ifndef FLOWBETWEENFLATPLATES_DIVERGENCE_H
#define FLOWBETWEENFLATPLATES_DIVERGENCE_H

#include "../../fields/scalarField/ScalarField.h"
#include "../../fields/vectorField/VectorField.h"


namespace fvc {

    ScalarField divergence(const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs);
};

#endif //FLOWBETWEENFLATPLATES_DIVERGENCE_H
