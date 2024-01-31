//
// Created by ruben on 31/01/24.
//

#ifndef SOLVERNS_WRITETURBULENTCHANNELFLOWDATA_H
#define SOLVERNS_WRITETURBULENTCHANNELFLOWDATA_H

#include <string>


void writeTurbulentChannelFlowData(double delta, double Lx, double Ly, double Lz, double Nx, double Ny, double Nz,
                                   double sx, double sy, double sz, double nu, double yPlusMin, double solverTolerance,
                                   std::string solver, int maxIter, double tInit, double tFinal, double f,
                                   double steadyStateCriterion, int writeInterval);

#endif //SOLVERNS_WRITETURBULENTCHANNELFLOWDATA_H
