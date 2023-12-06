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
    result.initialize(field1.field.size());

    for (int i = 0; i < field1.field.size(); ++i) {
        result.field[i] = field1.field[i] - field2.field[i];
    }
    return result;
}

// VectorField difference method
VectorField operator-(const VectorField &field1, const VectorField &field2) {
    VectorField result;
    result.initialize(field1.field.size());

    for (int i = 0; i < field1.field.size(); ++i) {
        result.field[i] = field1.field[i] - field2.field[i];
    }
    return result;
}

GeometricVector VectorField::operator[](int index) const {

    return field[index];
}