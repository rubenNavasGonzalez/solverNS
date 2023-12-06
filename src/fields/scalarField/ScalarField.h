//
// Created by ruben on 26/11/23.
//

#ifndef FLOWBETWEENFLATPLATES_SCALARFIELD_H
#define FLOWBETWEENFLATPLATES_SCALARFIELD_H

#include <vector>


class ScalarField {
public:
    // Scalar field
    std::vector<double> field;


    // ScalarField constructor and destructor
    ScalarField();
    ~ScalarField();


    // ScalarField methods
    void initialize(int length);
    double max();
    friend ScalarField operator-(const ScalarField& field1, const ScalarField& field2);
};


#endif //FLOWBETWEENFLATPLATES_SCALARFIELD_H
