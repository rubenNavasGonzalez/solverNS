//
// Created by ruben on 6/12/23.
//

#include "FvScalarEquation.h"
#include "../../math/linearSolver/BiCGSTAB/solveBiCGSTAB.h"


ScalarField FvScalarEquation::solve(const LinearSolverConfig& theLinearSolverConfig, const ScalarField& PhiOld) {

    if (theLinearSolverConfig.solver == "BiCGSTAB") {

        return solveBiCGSTAB(A, b, theLinearSolverConfig, PhiOld);
    } else {

        printf("ERROR. Incorrect solver selected!! \n");
        return PhiOld;
    }
}