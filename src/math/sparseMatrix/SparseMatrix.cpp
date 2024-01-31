//
// Created by ruben on 6/12/23.
//

#include "SparseMatrix.h"


// Default constructor and destructor
SparseMatrix::SparseMatrix() = default;
SparseMatrix::~SparseMatrix() = default;


// Matrix times vector operator
ScalarField operator*(const SparseMatrix &A, const ScalarField &field) {

    ScalarField result;
    result.resize( field.size() );

    if( !A.upperIndex.empty() ) {

        for(int i = 0; i < A.upperIndex.size(); i++)
        {

            result[A.upperIndex[i][0]] += A.upperValue[i]*field[A.upperIndex[i][1]];
        }
    }
    if( !A.diagIndex.empty() ) {

        for(int i = 0; i < A.diagIndex.size(); i++)
        {

            result[A.diagIndex[i][0]] += A.diagValue[i]*field[A.diagIndex[i][1]];
        }
    }
    if( !A.lowerIndex.empty() ) {

        for (int i = 0; i < A.lowerIndex.size(); i++) {

            result[A.lowerIndex[i][0]] += A.lowerValue[i]*field[A.lowerIndex[i][1]];
        }
    }


    return result;
}


// Sparse matrix multiplication (it is assumed A is a diagonal matrix)
SparseMatrix operator*(const SparseMatrix& A, const SparseMatrix& B) {

    // Initialize the product matrix
    SparseMatrix C;


    // Loop for all the B matrix diagonal elements
    for (int i = 0; i < B.diagIndex.size(); ++i) {

        C.diagValue.push_back( A.diagValue[i]*B.diagValue[i] );
        C.diagIndex.push_back({i,i});
    }


    // Loop for all the B matrix upper elements
    for (int i = 0; i < B.upperIndex.size(); ++i) {

        C.upperValue.push_back( A.diagValue[ B.upperIndex[i][0] ]*B.upperValue[i] );
        C.upperIndex.push_back( B.upperIndex[i] );
    }


    // Loop for all the B matrix lower elements
    for (int i = 0; i < B.lowerIndex.size(); ++i) {

        C.lowerValue.push_back( A.diagValue[ B.lowerIndex[i][0] ]*B.lowerValue[i] );
        C.lowerIndex.push_back( B.lowerIndex[i] );
    }


    return C;
}