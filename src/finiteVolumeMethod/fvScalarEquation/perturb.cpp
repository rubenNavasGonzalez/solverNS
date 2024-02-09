//
// Created by ruben on 7/02/24.
//

#include "FvScalarEquation.h"


void FvScalarEquation::perturb() {

    // Perturb the matrix to make the system of equation determined
    this->A.diagValue[0] = A.diagValue[0] * 1.1;
}