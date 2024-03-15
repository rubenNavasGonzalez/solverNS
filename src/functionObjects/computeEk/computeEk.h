//
// Created by ruben on 6/03/24.
//

#ifndef SOLVERNS_COMPUTEEK_H
#define SOLVERNS_COMPUTEEK_H

#include "../../fields/vectorField/VectorField.h"


double computeEk(const VectorField& u, const PolyMesh& theMesh, double t);

#endif //SOLVERNS_COMPUTEEK_H