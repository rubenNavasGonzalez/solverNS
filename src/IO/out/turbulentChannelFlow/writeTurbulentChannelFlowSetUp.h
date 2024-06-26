//
// Created by ruben on 31/01/24.
//

#ifndef SOLVERNS_WRITETURBULENTCHANNELFLOWSETUP_H
#define SOLVERNS_WRITETURBULENTCHANNELFLOWSETUP_H

#include <string>
#include "../../../turbulenceModeling/LES/modelLES.h"


void writeTurbulentChannelFlowSetUp(double delta, double Lx, double Ly, double Lz, double Nx, double Ny, double Nz,
                                   double sx, double sy, double sz, double nu, double yPlusMin, double solverTolerance,
                                   std::string solver, int maxIter, double tInit, double tFinal, double f, double steadyStateCriterion,
                                   double writeIntervalCSV,  double writeIntervalVTK, modelLES turbulenceModel, std::string filename);

#endif //SOLVERNS_WRITETURBULENTCHANNELFLOWSETUP_H
