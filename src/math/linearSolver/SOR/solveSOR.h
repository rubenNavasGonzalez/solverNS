//
// Created by ruben on 11/12/23.
//

#ifndef FLOWBETWEENFLATPLATES_SOLVESOR_H
#define FLOWBETWEENFLATPLATES_SOLVESOR_H

#include <string>
#include "../../../finiteVolumeMethod/fvScalarEquation/FvScalarEquation.h"


ScalarField solveSOR(const SparseMatrix& A, const ScalarField& b, const LinearSolverConfig& theLinearSolverConfig, const ScalarField& PhiOld);

#endif //FLOWBETWEENFLATPLATES_SOLVESOR_H
