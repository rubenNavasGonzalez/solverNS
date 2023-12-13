//
// Created by ruben on 6/12/23.
//

#include "solveBiCGSTAB.h"


ScalarField solveBiCGSTAB(const SparseMatrix& A, const ScalarField& b, const LinearSolverConfig& theLinearSolverConfig, const ScalarField& PhiOld) {

    int numberOfIterations = 0;

    ScalarField x = PhiOld;
    ScalarField x0 = PhiOld;
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

    if (r.maxAbs() <= theLinearSolverConfig.tolerance) {
        iterate = false;
        printf("\tSolution converged after %i iterations.", numberOfIterations);
    }

    while (iterate) {
        rho = r0*rPrev;

        if(numberOfIterations != 0)
        {
            beta = (rho/rhoPrev)*(alpha/omegaPrev);
            p = rPrev + beta*( pPrev - omegaPrev*vPrev );
        } else {
            p = rPrev;
        }

        v = A*p;
        alpha = rho / (r0*v);
        s = rPrev - alpha*v;
        t = A*s;
        omega = (t * s)/(t * t);

        x = x0 + alpha*p + omega*s;
        r = b - A*x;
        numberOfIterations++;

        if(r.maxAbs() <= theLinearSolverConfig.tolerance) {
            iterate = false;
            printf("\tSolution converged after %i iterations.", numberOfIterations);
        } else if(numberOfIterations > theLinearSolverConfig.maxIter) {
            iterate = false;
            printf("\tMaximum number of iterations reached.");
        } else {
            r = s - omega*t;
            rhoPrev = rho;
            omegaPrev = omega;
            x0 = x;
            rPrev = r;
            pPrev = p;
            vPrev = v;
            sPrev = s;
        }
    }

    return x;
}