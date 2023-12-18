//
// Created by ruben on 18/12/23.
//

#ifndef FLOWBETWEENFLATPLATES_INTERPOLATEVECTORFIELDFROMELEMENTS2NODES_H
#define FLOWBETWEENFLATPLATES_INTERPOLATEVECTORFIELDFROMELEMENTS2NODES_H

#include "../../fields/vectorField/VectorField.h"
#include "../../boundaryConditions/vectorBoundaryConditions/VectorBoundaryConditions.h"


VectorField interpolateVectorFieldFromElements2Nodes(const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs);

#endif //FLOWBETWEENFLATPLATES_INTERPOLATEVECTORFIELDFROMELEMENTS2NODES_H
