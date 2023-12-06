//
// Created by ruben on 6/12/23.
//

#include "FvScalarEquation.h"


// FvScalarEquation default constructor and destructor
FvScalarEquation::FvScalarEquation() = default;
FvScalarEquation::~FvScalarEquation() = default;


// FvScalarEquation overload == operator
FvScalarEquation operator==(const SparseMatrix& A, const ScalarField& b) {

    FvScalarEquation result;

    result.A = A;
    result.b = b;


    return result;
}
