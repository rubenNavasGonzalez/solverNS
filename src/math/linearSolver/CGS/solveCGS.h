//
// Created by ruben on 26/01/24.
//

#ifndef SOLVERNS_SOLVECGS_H
#define SOLVERNS_SOLVECGS_H

#include <string>
#include "../../../finiteVolumeMethod/fvScalarEquation/FvScalarEquation.h"


ScalarField solveCGS(const SparseMatrix& A, const ScalarField& b, const LinearSolverConfig& theLinearSolverConfig, const ScalarField& PhiOld);

#endif //SOLVERNS_SOLVECGS_H
