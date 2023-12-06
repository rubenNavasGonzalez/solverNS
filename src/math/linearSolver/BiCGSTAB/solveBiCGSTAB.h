//
// Created by ruben on 6/12/23.
//

#ifndef FLOWBETWEENFLATPLATES_SOLVEBICGSTAB_H
#define FLOWBETWEENFLATPLATES_SOLVEBICGSTAB_H

#include <string>
#include "../../../finiteVolumeMethod/fvScalarEquation/FvScalarEquation.h"


ScalarField solveBiCGSTAB(const SparseMatrix& A, const ScalarField& b, const LinearSolverConfig& theLinearSolverConfig, const ScalarField& PhiOld);

#endif //FLOWBETWEENFLATPLATES_SOLVEBICGSTAB_H
