//
// Created by ruben on 10/04/24.
//

#ifndef SOLVERNS_INITIALIZEFROMDIFFERENTMESH_H
#define SOLVERNS_INITIALIZEFROMDIFFERENTMESH_H

#include "ScalarField.h"


ScalarField interpolateScalarFieldFromDifferentMesh(const PolyMesh& theMesh, double t);

#endif //SOLVERNS_INITIALIZEFROMDIFFERENTMESH_H
