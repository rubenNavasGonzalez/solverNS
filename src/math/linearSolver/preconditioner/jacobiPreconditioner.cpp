//
// Created by ruben on 31/01/24.
//

#include "jacobiPreconditioner.h"


SparseMatrix jacobiPreconditioner(const SparseMatrix& A) {

    // Initialize the preconditioner matrix
    SparseMatrix P = A;


    // Assemble the preconditioner matrix
    for (int i = 0; i < P.diagIndex.size(); ++i) {

        P.diagValue[i] = 1/P.diagValue[i];
    }


    return P;
}