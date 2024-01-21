//
// Created by ruben on 28/12/23.
//

#ifndef SOLVERNS_WRITEPRESSUREVELOCITY2TXT_H
#define SOLVERNS_WRITEPRESSUREVELOCITY2TXT_H

#include "../fields/scalarField/ScalarField.h"
#include "../fields/vectorField/VectorField.h"


void writePressureVelocity2TXT(const PolyMesh& theMesh, const ScalarField& p, const VectorField& u, double t);

#endif //SOLVERNS_WRITEPRESSUREVELOCITY2TXT_H
