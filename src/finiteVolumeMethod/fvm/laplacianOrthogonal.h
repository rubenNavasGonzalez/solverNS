//
// Created by ruben on 6/12/23.
//

//#ifndef FLOWBETWEENFLATPLATES_LAPLACIANORTHOGONAL_H
//#define FLOWBETWEENFLATPLATES_LAPLACIANORTHOGONAL_H

#include "../../math/sparseMatrix/SparseMatrix.h"
#include "../../mesh/polyMesh/PolyMesh.h"
#include "../../boundaryConditions/scalarBoundaryConditions/ScalarBoundaryConditions.h"


namespace fvm {

    SparseMatrix laplacianOrthogonal(const PolyMesh& theMesh);
}

//#endif //FLOWBETWEENFLATPLATES_LAPLACIANORTHOGONAL_H
