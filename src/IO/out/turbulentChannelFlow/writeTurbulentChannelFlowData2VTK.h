//
// Created by ruben on 22/01/24.
//

#ifndef SOLVERNS_WRITETURBULENTCHANNELFLOWDATA2VTK_H
#define SOLVERNS_WRITETURBULENTCHANNELFLOWDATA2VTK_H

#include "../../../fields/scalarField/ScalarField.h"
#include "../../../fields/vectorField/VectorField.h"


void writeTurbulentChannelFlowData2VTK(const PolyMesh& theMesh, const ScalarField& p, const VectorField& u, const VectorField& omega,
                                       const ScalarField& nut, const ScalarField& QCrit, const ScalarBoundaryConditions& pBCs,
                                       const VectorBoundaryConditions& uBCs, const VectorBoundaryConditions& omegaBCs,
                                       const ScalarBoundaryConditions& nutBCs, const ScalarBoundaryConditions& QCritBCs, double t);

#endif //SOLVERNS_WRITETURBULENTCHANNELFLOWDATA2VTK_H