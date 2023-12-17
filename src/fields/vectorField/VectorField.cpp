//
// Created by ruben on 26/11/23.
//

#include "VectorField.h"


// VectorField default constructor and destructor
VectorField::VectorField() = default;
VectorField::~VectorField() = default;


// VectorField sum operator
VectorField operator+(const VectorField &field1, const VectorField &field2) {

    VectorField result;
    result.assign(field1.size(), {0,0,0});

    for (int i = 0; i < field1.size(); ++i) {

        result[i] = field1[i] + field2[i];
    }


    return result;
}


// VectorField difference method
VectorField operator-(const VectorField &field1, const VectorField &field2) {

    VectorField result;
    result.assign(field1.size(), {0,0,0});

    for (int i = 0; i < field1.size(); ++i) {
        result[i] = field1[i] - field2[i];
    }


    return result;
}

VectorField operator*(const double &k, const VectorField &field) {

    VectorField result = field;

    for (int i = 0; i < field.size(); ++i) {

        result[i] = k*result[i];
    }


    return result;
}