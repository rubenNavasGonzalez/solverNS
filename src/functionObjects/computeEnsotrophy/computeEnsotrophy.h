//
// Created by ruben on 6/03/24.
//

#ifndef SOLVERNS_COMPUTEENSOTROPHY_H
#define SOLVERNS_COMPUTEENSOTROPHY_H


#include "../../fields/vectorField/VectorField.h"


double computeEnsotrophy(const VectorField& omega, const PolyMesh& theMesh, double t);

#endif //SOLVERNS_COMPUTEENSOTROPHY_H
