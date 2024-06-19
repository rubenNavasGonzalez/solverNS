//
// Created by ruben on 10/04/24.
//

#ifndef SOLVERNS_WRITEDIFFERENTIALLYHEATEDCAVITYSETUP_H
#define SOLVERNS_WRITEDIFFERENTIALLYHEATEDCAVITYSETUP_H

#include <string>
#include "../../../turbulenceModeling/LES/modelLES.h"


void writeDifferentiallyHeatedCavitySetUp(double LRef, double Lx, double Ly, double Lz, double Nx, double Ny, double Nz,
                                          double sx, double sy, double sz, double deltaXMin, double deltaYMin, double Pr, double Ra,
                                          double solverTolerance, std::string solver, int maxIter, double tInit, double tFinal,
                                          double DeltaT, double steadyStateCriterion, double writeIntervalCSV,  double writeIntervalVTK,
                                          modelLES turbulenceModel, std::string filename);

#endif //SOLVERNS_WRITEDIFFERENTIALLYHEATEDCAVITYSETUP_H
