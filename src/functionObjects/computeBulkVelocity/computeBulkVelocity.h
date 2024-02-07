//
// Created by ruben on 26/12/23.
//

#ifndef SOLVERNS_COMPUTEBULKVELOCITY_H
#define SOLVERNS_COMPUTEBULKVELOCITY_H

#include "../../fields/vectorField/VectorField.h"


double computeBulkVelocity(const VectorField& u, const PolyMesh& theMesh, const VectorBoundaryConditions& uBCs, int k, double t);

#endif //SOLVERNS_COMPUTEBULKVELOCITY_H
