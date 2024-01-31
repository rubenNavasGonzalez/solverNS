//
// Created by ruben on 6/12/23.
//

#include "FvScalarEquation.h"
#include "../../math/linearSolver/CGS/solveCGS.h"
#include "../../math/linearSolver/BiCGSTAB/solveBiCGSTAB.h"
#include "../../math/linearSolver/SOR/solveSOR.h"



ScalarField FvScalarEquation::solve(const LinearSolverConfig& theLinearSolverConfig, const ScalarField& PhiOld) {

    if (theLinearSolverConfig.solver == "CGS") {

        return solveCGS(A, b, theLinearSolverConfig, PhiOld);
    } else if (theLinearSolverConfig.solver == "BiCGSTAB") {

        return solveBiCGSTAB(A, b, theLinearSolverConfig, PhiOld);
    } else if (theLinearSolverConfig.solver == "SOR") {

        return solveSOR(A, b, theLinearSolverConfig, PhiOld);
    } else {

        printf("ERROR. Incorrect solver selected!! \n");
        return PhiOld;
    }
}