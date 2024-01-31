//
// Created by ruben on 6/12/23.
//

#ifndef FLOWBETWEENFLATPLATES_SPARSEMATRIX_H
#define FLOWBETWEENFLATPLATES_SPARSEMATRIX_H

#include <vector>
#include "../../fields/scalarField/ScalarField.h"


class SparseMatrix {
public:
    // SparseMatrix attributes
    int nRows{}, nCols{};
    std::vector<double> upperValue, diagValue, lowerValue;
    std::vector< std::vector<int> > upperIndex, diagIndex, lowerIndex;


    // SparseMatrix constructor and destructor
    SparseMatrix();
    ~SparseMatrix();


    // SparseMatrix methods
    std::vector<double> getValue(std::vector<int> ij);
    void addValue(double value, std::vector<int> ij);
    friend ScalarField operator*(const SparseMatrix& A, const ScalarField& field);
    friend SparseMatrix operator*(const SparseMatrix& A, const SparseMatrix& B);
};


#endif //FLOWBETWEENFLATPLATES_SPARSEMATRIX_H
