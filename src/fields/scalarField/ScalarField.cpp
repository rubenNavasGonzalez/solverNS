//
// Created by ruben on 26/11/23.
//

#include "ScalarField.h"


// ScalarField default constructor and destructor
ScalarField::ScalarField() = default;
ScalarField::~ScalarField() = default;


// Difference of ScalarField
ScalarField operator-(const ScalarField &field1, const ScalarField &field2) {
    ScalarField result;

    for (int i = 0; i < field1.field.size(); ++i) {
        result.field.push_back( field1.field[i] - field2.field[i] );
    }
    return result;
}


ScalarField operator+(const ScalarField &field1, const ScalarField &field2) {
    ScalarField result;

    for (int i = 0; i < field1.field.size(); ++i) {
        result.field.push_back( field1[i] + field2[i] );
    }
    return result;
}


// Dot product operator of ScalarField
double operator*(const ScalarField &field1, const ScalarField &field2) {
    double result = 0;

    for (int i = 0; i < field1.field.size(); ++i) {
        result += field1[i]*field2[i];
    }

    return result;
}


ScalarField operator*(const double &k, const ScalarField &field) {
    ScalarField result;

    for (int i = 0; i < field.field.size(); ++i) {
        result.field.push_back( k*field[i] );
    }
    return result;
}


// Operator [] overload
double ScalarField::operator[](int index) const {

    return field[index];
}