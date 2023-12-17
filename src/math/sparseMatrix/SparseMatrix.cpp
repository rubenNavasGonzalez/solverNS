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