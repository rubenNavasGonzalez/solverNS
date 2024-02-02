//
// Created by ruben on 26/01/24.
//

#include "solveCGS.h"


ScalarField solveCGS(const SparseMatrix& A, const ScalarField& b, const LinearSolverConfig& theLinearSolverConfig, const ScalarField& PhiOld) {

    // Auxiliary variables initialization
    int numberOfIterations = 0;
    ScalarField r, rOld, rHat, p, t, q, u, Phi;
    double beta, alpha;
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
    u = r;

    // Iterative process until reach convergence or exceed the maximum number of iterations
    while (iterate) {


        // Update the solution
        t = A * p;
        alpha = (r * rHat) / (rHat * t);
        q = u - alpha * t;
        Phi = Phi + alpha * (u + q);


        // Update the residual
        r = rOld - alpha * (A * (u + q));


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
            beta = (r * rHat) / (rOld * rHat);
            u = r + beta * q;
            p = u + beta * (q + beta * p);

            rOld = r;
        }
    }


    return Phi;
}