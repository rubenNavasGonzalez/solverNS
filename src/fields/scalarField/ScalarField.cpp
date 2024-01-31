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
    result.assign(field1.size(), 0);

    for (int i = 0; i < field1.size(); ++i) {

        result[i] = field1[i] - field2[i] ;
    }


    return result;
}


ScalarField operator+(const ScalarField &field1, const ScalarField &field2) {

    ScalarField result;
    result.assign(field1.size(), 0);

    for (int i = 0; i < field1.size(); ++i) {

        result[i] = field1[i] + field2[i] ;
    }


    return result;
}


// Dot product operator of ScalarField
double operator*(const ScalarField &field1, const ScalarField &field2) {

    double result = 0;

    for (int i = 0; i < field1.size(); ++i) {

        result += field1[i]*field2[i];
    }


    return result;
}


ScalarField operator*(const double &k, const ScalarField &field) {

    ScalarField result;
    result.assign(field.size(), 0);


    for (int i = 0; i < field.size(); ++i) {

        result[i] = k*field[i];
    }


    return result;
}


ScalarField operator/(const ScalarField& field, const double& k) {

    ScalarField result = field;

    for (int i = 0; i < field.size(); ++i) {

        result[i] = result[i]/k;
    }


    return result;
}
