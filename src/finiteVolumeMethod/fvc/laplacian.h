//
// Created by ruben on 26/11/23.
//

#ifndef FLOWBETWEENFLATPLATES_LAPLACIAN_H
#define FLOWBETWEENFLATPLATES_LAPLACIAN_H

#include "divergence.h"


namespace fvc {

    VectorField laplacian(const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs);
}


#endif //FLOWBETWEENFLATPLATES_LAPLACIAN_H
