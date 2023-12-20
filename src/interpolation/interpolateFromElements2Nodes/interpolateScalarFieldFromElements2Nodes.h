//
// Created by ruben on 18/12/23.
//

#ifndef FLOWBETWEENFLATPLATES_INTERPOLATESCALARFIELDFROMELEMENTS2NODES_H
#define FLOWBETWEENFLATPLATES_INTERPOLATESCALARFIELDFROMELEMENTS2NODES_H

#include "../../fields/scalarField/ScalarField.h"
#include "../../mesh/polyMesh/PolyMesh.h"
#include "../../boundaryConditions/scalarBoundaryConditions/ScalarBoundaryConditions.h"


ScalarField interpolateScalarFieldFromElements2Nodes(const ScalarField& Phi, const PolyMesh& theMesh, const ScalarBoundaryConditions& PhiBCs);

#endif //FLOWBETWEENFLATPLATES_INTERPOLATESCALARFIELDFROMELEMENTS2NODES_H
