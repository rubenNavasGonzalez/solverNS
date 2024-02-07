//
// Created by ruben on 26/11/23.
//

#include "ScalarField.h"


// ScalarField default constructor and destructor
ScalarField::ScalarField() = default;
ScalarField::~ScalarField() = default;



// ScalarField difference
ScalarField operator-(const ScalarField &field1, const ScalarField &field2) {

    // Initialize the result
    ScalarField result;
    result.assign(field1.size(), 0);


    // Perform the operation
    for (int i = 0; i < field1.size(); ++i) {

        result[i] = field1[i] - field2[i] ;
    }


    return result;
}



// ScalarField sum
ScalarField operator+(const ScalarField &field1, const ScalarField &field2) {

    // Initialize the result
    ScalarField result;
    result.assign(field1.size(), 0);


    // Perform the operation
    for (int i = 0; i < field1.size(); ++i) {

        result[i] = field1[i] + field2[i] ;
    }


    return result;
}



// ScalarField dot product
double operator*(const ScalarField &field1, const ScalarField &field2) {

    // Initialize the result
    double result = 0;


    // Perform the operation
    for (int i = 0; i < field1.size(); ++i) {

        result += field1[i] * field2[i];
    }


    return result;
}



// Constant time ScalarField
ScalarField operator*(const double &k, const ScalarField &field) {

    // Initialize the result
    ScalarField result;
    result.assign(field.size(), 0);


    // Perform the operation
    for (int i = 0; i < field.size(); ++i) {

        result[i] = k * field[i];
    }


    return result;
}



// ScalarField divided by a constant
ScalarField operator/(const ScalarField& field, const double& k) {

    // Initialize the result
    ScalarField result;
    result.assign(field.size(), 0);


    // Perform the operation
    for (int i = 0; i < field.size(); ++i) {

        result[i] = field[i] / k;
    }


    return result;
}