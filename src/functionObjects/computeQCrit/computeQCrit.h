//
// Created by ruben on 31/01/24.
//

#ifndef SOLVERNS_COMPUTEQCRIT_H
#define SOLVERNS_COMPUTEQCRIT_H

#include "../../fields/scalarField/ScalarField.h"
#include "../../fields/tensorField/TensorField.h"


ScalarField computeQCrit(const PolyMesh& theMesh, TensorField Phi);

#endif //SOLVERNS_COMPUTEQCRIT_H
