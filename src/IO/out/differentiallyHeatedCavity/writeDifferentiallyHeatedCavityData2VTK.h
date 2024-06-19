//
// Created by ruben on 6/04/24.
//

#ifndef SOLVERNS_WRITEDIFFERENTIALLYHEATEDCAVITYDATA2VTK_H
#define SOLVERNS_WRITEDIFFERENTIALLYHEATEDCAVITYDATA2VTK_H

#include "../../../fields/scalarField/ScalarField.h"
#include "../../../fields/vectorField/VectorField.h"


void writeDifferentiallyHeatedCavityData2VTK(const PolyMesh& theMesh, const ScalarField& p, const ScalarField& T, const VectorField& u,
                                       const VectorField& omega, const ScalarField& nut, const ScalarField& QCrit,
                                       const ScalarBoundaryConditions& pBCs, const ScalarBoundaryConditions& TBCs,
                                       const VectorBoundaryConditions& uBCs, const VectorBoundaryConditions& omegaBCs,
                                       const ScalarBoundaryConditions& nutBCs, const ScalarBoundaryConditions& QCritBCs, double t);

#endif //SOLVERNS_WRITEDIFFERENTIALLYHEATEDCAVITYDATA2VTK_H