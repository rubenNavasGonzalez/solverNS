//
// Created by ruben on 11/12/23.
//

#include <iostream>
#include "solveSOR.h"


ScalarField solveSOR(const SparseMatrix& A, const ScalarField& b, const LinearSolverConfig& theLinearSolverConfig, const ScalarField& PhiOld) {

    double maxResidual;
    int numberOfIterations = 0;

    ScalarField x0 = PhiOld;
    ScalarField x = PhiOld;

    double xi, aii;

    bool iterate = true;
    while (iterate) {
        for (int i = 0; i < x.field.size(); ++i) {
            xi = 0;
            if ( !A.upperValue.empty() ) {
                for (int j = 0; j < A.upperValue.size(); ++j) {
                    if ( A.upperIndex[j][0] == i ) {
                        xi += A.upperValue[j]*x.field[ A.upperIndex[j][1] ];
                    }
                }
            }
            if ( !A.diagValue.empty() ) {
                for (int j = 0; j < A.diagValue.size(); ++j) {
                    if ( A.diagIndex[j][0] == i ) {
                        aii = A.diagValue[j];
                    }
                }
            }
            if ( !A.lowerValue.empty() ) {
                for (int j = 0; j < A.lowerValue.size(); ++j) {
                    if ( A.lowerIndex[j][0] == i ) {
                        xi += A.lowerValue[j]*x.field[ A.lowerIndex[j][1] ];
                    }
                }
            }

            xi = 1/aii*(b.field[i] - xi);
            x.field[i] = xi;
        }
        numberOfIterations++;

        ScalarField r = A*x - b;

        if(r.maxAbs() <= theLinearSolverConfig.tolerance)
        {
            iterate = false;
            std::cout << "\tSolution converged after " << numberOfIterations << " iterations.\n";
        }
        else if(numberOfIterations > theLinearSolverConfig.maxIter)
        {
            iterate = false;
            std::cout << "\tMaximum number of iterations reached.\n";
        }
        else
        {
            x0 = x;
        }
    }

    return x;
}