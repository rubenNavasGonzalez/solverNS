//
// Created by ruben on 27/11/23.
//

#ifndef FLOWBETWEENFLATPLATES_LAPLACIANORTHOGONAL_H
#define FLOWBETWEENFLATPLATES_LAPLACIANORTHOGONAL_H

#include "../../fields/vectorField/VectorField.h"
#include "../../mesh/polyMesh/PolyMesh.h"
#include "../../boundaryConditions/vectorBoundaryConditions/VectorBoundaryConditions.h"


namespace fvc {

    VectorField laplacianOrthogonal(const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs);
}

#endif //FLOWBETWEENFLATPLATES_LAPLACIANORTHOGONAL_H
