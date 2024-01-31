//
// Created by ruben on 26/01/24.
//

#include "solveCGS.h"
#include "../preconditioner/jacobiPreconditioner.h"


ScalarField solveCGS(const SparseMatrix& A, const ScalarField& b, const LinearSolverConfig& theLinearSolverConfig, const ScalarField& PhiOld) {

    // Auxiliary variables initialization
    int numberOfIterations = 0;
    ScalarField r, rOld, p, pOld, Phi, PhiAux;
    double beta, alpha;
    bool iterate = true;


    // Check if the Phi field is already converged
    PhiAux = PhiOld;
    r = b - A*PhiAux;
    rOld = r;

    if (r.maxAbs() <= theLinearSolverConfig.tolerance) {

        printf("\tSolution converged after %i iterations. \n", numberOfIterations);
        return PhiAux;
    }


    // Compute the initial search direction
    p = r;
    pOld = p;

    // Iterative process until reach convergence or exceed the maximum number of iterations
    while (iterate) {


        // Update the solution
        alpha = (rOld*rOld) / (pOld*(A*pOld));
        Phi = PhiAux + alpha*pOld;


        // Update the residual
        r = rOld - alpha*(A*pOld);


        // Compute the search direction
        beta = (r*r) / (rOld*rOld);
        p = r + beta*pOld;


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

            rOld = r;
            pOld = p;
            PhiAux = Phi;
        }
    }


    return Phi;
}