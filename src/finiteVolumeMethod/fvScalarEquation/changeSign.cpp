//
// Created by ruben on 9/02/24.
//

#include "FvScalarEquation.h"


void FvScalarEquation::changeSign() {

    // Change the sign of the matrix of the system of equations
    for (int i = 0; i < A.diagValue.size(); ++i) {

        this->A.diagValue[i] = A.diagValue[i] * (-1);
    }

    for (int i = 0; i < A.lowerValue.size(); ++i) {

        this->A.lowerValue[i] = A.lowerValue[i] * (-1);
    }

    for (int i = 0; i < this->A.lowerValue.size(); ++i) {

        this->A.upperValue[i] = A.upperValue[i] * (-1);
    }


    // Change the sign of the vector of the system of equations
    for (int i = 0; i < b.size(); ++i) {

        b[i] = b[i] * (-1);
    }
}