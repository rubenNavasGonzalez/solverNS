//
// Created by ruben on 1/02/24.
//

#ifndef SOLVERNS_SOLVECG_H
#define SOLVERNS_SOLVECG_H

#include "../../../finiteVolumeMethod/fvScalarEquation/FvScalarEquation.h"


ScalarField solveCG(const SparseMatrix& A, const ScalarField& b, const LinearSolverConfig& theLinearSolverConfig, const ScalarField& PhiOld);

#endif //SOLVERNS_SOLVECG_H