//
// Created by ruben on 28/12/23.
//

#ifndef SOLVERNS_WRITEPRESSUREVELOCITY2VTK_H
#define SOLVERNS_WRITEPRESSUREVELOCITY2VTK_H

#include "../fields/scalarField/ScalarField.h"
#include "../fields/vectorField/VectorField.h"


void writePressureVelocity2VTK(const PolyMesh& theMesh, const ScalarField& p, const VectorField& u, const ScalarBoundaryConditions& pBCs,
                               const VectorBoundaryConditions& uBCs, double t);

#endif //SOLVERNS_WRITEPRESSUREVELOCITY2VTK_H
