//
// Created by ruben on 6/12/23.
//

#include "SparseMatrix.h"


// Default constructor and destructor
SparseMatrix::SparseMatrix() = default;
SparseMatrix::~SparseMatrix() = default;



// Matrix times vector operator (it is assumed A is a symmetric matrix)
ScalarField operator*(const SparseMatrix &A, const ScalarField &field) {

    // Initialize the result
    ScalarField result;
    result.assign(field.size(), 0);


    // Computations regarding the upper and lower triangular part of A
    for (int i = 0; i < A.upperIndex.size(); i++) {

        // Upper triangular part of A
        result[A.upperIndex[i][0]] += A.upperValue[i] * field[A.upperIndex[i][1]];

        // Lower triangular part of A
        result[A.upperIndex[i][1]] += A.upperValue[i] * field[A.upperIndex[i][0]];
    }


    // Computations regarding the diagonal part of A
    for (int i = 0; i < A.diagIndex.size(); i++) {

        result[A.diagIndex[i][0]] += A.diagValue[i] * field[A.diagIndex[i][1]];
    }


    return result;
}



// Sparse matrix multiplication (it is assumed A is a diagonal matrix and B a symmetric matrix)
SparseMatrix operator*(const SparseMatrix& A, const SparseMatrix& B) {

    // Initialize result
    SparseMatrix C;


    // Computations regarding the diagonal part of A
    for (int i = 0; i < B.diagIndex.size(); ++i) {

        C.diagValue.push_back( A.diagValue[0] * B.diagValue[i] );
        C.diagIndex.push_back( {i, i} );
    }


    // Computations regarding the upper and lower triangular part of A
    for (int i = 0; i < B.upperIndex.size(); ++i) {

        // Upper triangular part of A
        C.upperValue.push_back( A.diagValue[0] * B.upperValue[i] );
        C.upperIndex.push_back( B.upperIndex[i] );

        // Lower triangular part of A
        C.lowerValue.push_back( A.diagValue[0] * B.upperValue[i] );
        C.lowerIndex.push_back( {B.upperIndex[i][1], B.upperIndex[i][0]} );
    }

    return C;
}