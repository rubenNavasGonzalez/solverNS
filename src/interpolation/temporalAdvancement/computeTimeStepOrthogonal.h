//
// Created by ruben on 15/12/23.
//

#ifndef FLOWBETWEENFLATPLATES_COMPUTETIMESTEPORTHOGONAL_H
#define FLOWBETWEENFLATPLATES_COMPUTETIMESTEPORTHOGONAL_H

#include "../../mesh/polyMesh/PolyMesh.h"
#include "../../fields/vectorField/VectorField.h"


double computeTimeStepOrthogonal(const PolyMesh& theMesh, const VectorField& u, double nu, double f);

#endif //FLOWBETWEENFLATPLATES_COMPUTETIMESTEPORTHOGONAL_H
