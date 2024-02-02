//
// Created by ruben on 1/02/24.
//

#include "solveCG.h"
#include "../preconditioner/jacobiPreconditioner.h"


ScalarField solveCG(const SparseMatrix& A, const ScalarField& b, const LinearSolverConfig& theLinearSolverConfig, const ScalarField& PhiOld) {

    // Auxiliary variables initialization
    int numberOfIterations = 0;
    ScalarField r, rOld, p, t, Phi;
    double beta, alpha;
    bool iterate = true;


    // Check if the Phi field is already converged
    Phi = PhiOld;
    r = b - A*Phi;
    rOld = r;

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
        alpha = (r * r) / (p * t);
        Phi = Phi + alpha * p;


        // Update the residual
        r = r - alpha * t;


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
            beta = (r * r) / (rOld * rOld);
            p = r + beta * p;

            rOld = r;
        }
    }


    return Phi;
}
