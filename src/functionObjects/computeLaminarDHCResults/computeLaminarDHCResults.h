//
// Created by ruben on 6/04/24.
//

#ifndef SOLVERNS_COMPUTELAMINARDHCRESULTS_H
#define SOLVERNS_COMPUTELAMINARDHCRESULTS_H

#include "../../fields/vectorField/VectorField.h"
#include "../../fields/scalarField/ScalarField.h"


void computeLaminarDHCResults(const PolyMesh& theMesh, int Nx, int Ny, int Nz, double Lx, double Ly, double t, const VectorField& u,
                              const ScalarField& T, const ScalarBoundaryConditions& TBCs);

#endif //SOLVERNS_COMPUTELAMINARDHCRESULTS_H
