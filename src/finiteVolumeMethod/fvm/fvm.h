//
// Created by ruben on 7/12/23.
//

#ifndef FLOWBETWEENFLATPLATES_FVM_H
#define FLOWBETWEENFLATPLATES_FVM_H

#include "../../math/sparseMatrix/SparseMatrix.h"
#include "../../mesh/polyMesh/PolyMesh.h"
#include "../../boundaryConditions/scalarBoundaryConditions/ScalarBoundaryConditions.h"


namespace fvm {

    SparseMatrix laplacianOrthogonal(const PolyMesh& theMesh);
}

#endif //FLOWBETWEENFLATPLATES_FVM_H
