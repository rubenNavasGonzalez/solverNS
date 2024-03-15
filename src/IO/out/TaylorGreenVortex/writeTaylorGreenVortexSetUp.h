//
// Created by ruben on 3/03/24.
//

#ifndef SOLVERNS_WRITETAYLORGREENVORTEXSETUP_H
#define SOLVERNS_WRITETAYLORGREENVORTEXSETUP_H

#include <string>
#include "../../../turbulenceModeling/LES/modelLES.h"


void writeTaylorGreenVortexSetUp(double L, double Lx, double Ly, double Lz, int Nx, int Ny, int Nz, double sx, double sy,
                                 double sz, double nu, double solverTolerance, std::string solver, int maxIter, double tInit,
                                 double tFinal, double DeltaT, double steadyStateCriterion, double writeIntervalCSV,  double writeIntervalVTK,
                                 modelLES turbulenceModel, std::string filename);

#endif //SOLVERNS_WRITETAYLORGREENVORTEXSETUP_H
