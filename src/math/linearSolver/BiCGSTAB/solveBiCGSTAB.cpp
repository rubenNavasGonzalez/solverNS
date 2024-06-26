//
// Created by ruben on 6/12/23.
//

#include "solveBiCGSTAB.h"


ScalarField solveBiCGSTAB(const SparseMatrix& A, const ScalarField& b, const LinearSolverConfig& theLinearSolverConfig, const ScalarField& PhiOld) {


     // Auxiliary variables initialization
     int numberOfIterations = 0;
     ScalarField r, rOld, rHat, p, s, t, q, Phi;
     double alpha, omega, beta;
     bool iterate = true;


    // Check if the Phi field is already converged
    Phi = PhiOld;
    r = b - A*Phi;
    rOld = r;
    rHat = r;

    if (r.maxAbs() <= theLinearSolverConfig.tolerance) {

        printf("\tSolution converged after %i iterations. \n", numberOfIterations);

        return Phi;
    }


    // Compute the initial search direction
    p = r;


    // Iterative process until reach convergence or exceed the maximum number of iterations
    while (iterate) {


        // Update the solution
        t = A * p;
        alpha = (r * rHat) / (t * rHat);
        s = r - alpha * t;
        q = A * s;
        omega = (q * s) / (q * q);
        Phi = Phi + (alpha * p) + (omega * s);
        r = s - omega * q;


        // Update the number of iterations
        numberOfIterations++;


        // Check convergence
        if(r.maxAbs() <= theLinearSolverConfig.tolerance) {

            iterate = false;
            printf("\tSolution converged after %i iterations. \n", numberOfIterations);
        } else if(numberOfIterations > theLinearSolverConfig.maxIter) {

            iterate = false;
            printf("\tMaximum number of iterations reached. \n");
        } else {

            // Compute the search direction and update parameters
            beta = (alpha / omega) * (r * rHat) / (rOld * rHat);
            p = r + beta * (p - omega * t);

            rOld = r;
        }
    }


    return Phi;
}