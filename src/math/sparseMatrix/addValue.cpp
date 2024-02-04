//
// Created by ruben on 6/12/23.
//

#include "SparseMatrix.h"


void SparseMatrix::addValue(double value, std::vector<int> ij) {

    if(ij[0] < ij[1]) {

        upperValue.push_back(value);
        upperIndex.push_back(ij);
    } else if(ij[0] > ij[1]) {

        lowerValue.push_back(value);
        lowerIndex.push_back(ij);
    } else {
        diagValue.push_back(value);
        diagIndex.push_back(ij);
    }
}