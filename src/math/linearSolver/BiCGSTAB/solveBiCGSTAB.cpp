//
// Created by ruben on 6/12/23.
//

#include "solveBiCGSTAB.h"


ScalarField solveBiCGSTAB(const SparseMatrix& A, const ScalarField& b, const LinearSolverConfig& theLinearSolverConfig, const ScalarField& PhiOld) {

    double maxResidual;
    int numberOfIterations = 0;

    ScalarField x = PhiOld;
    ScalarField r;
    ScalarField rPrev;
    ScalarField r0;
    ScalarField p;
    ScalarField pPrev;
    ScalarField v;
    ScalarField vPrev;
    ScalarField s;
    ScalarField sPrev;
    ScalarField t;

    double rho;
    double rhoPrev;
    double beta;
    double alpha;
    double omega;
    double omegaPrev;

    bool iterate = true;

    r = b - A*x;
    rPrev = r;
    r0 = r;


}