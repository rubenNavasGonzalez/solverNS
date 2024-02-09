//
// Created by ruben on 21/01/24.
//

#include "Tensor.h"


// Tensor default constructor and destructor
Tensor::Tensor() {

    this->resize(3, {0,0,0});
}

Tensor::~Tensor() = default;


// Tensor sum operator
Tensor operator+(const Tensor &T1, const Tensor &T2) {

    Tensor result;

    result[0][0] = T1[0][0] + T2[0][0];
    result[0][1] = T1[0][1] + T2[0][1];
    result[0][2] = T1[0][2] + T2[0][2];
    result[1][0] = T1[1][0] + T2[1][0];
    result[1][1] = T1[1][1] + T2[1][1];
    result[1][2] = T1[1][2] + T2[1][2];
    result[2][0] = T1[2][0] + T2[2][0];
    result[2][1] = T1[2][1] + T2[2][1];
    result[2][2] = T1[2][2] + T2[2][2];


    return result;
}


// Tensor subtract operator
Tensor operator-(const Tensor &T1, const Tensor &T2) {

    Tensor result;

    result[0][0] = T1[0][0] - T2[0][0];
    result[0][1] = T1[0][1] - T2[0][1];
    result[0][2] = T1[0][2] - T2[0][2];
    result[1][0] = T1[1][0] - T2[1][0];
    result[1][1] = T1[1][1] - T2[1][1];
    result[1][2] = T1[1][2] - T2[1][2];
    result[2][0] = T1[2][0] - T2[2][0];
    result[2][1] = T1[2][1] - T2[2][1];
    result[2][2] = T1[2][2] - T2[2][2];


    return result;
}


// Tensor scalar times a tensor operator
Tensor operator*(const double &k, const Tensor &T) {

    Tensor result;

    result[0][0] = k*T[0][0];
    result[0][1] = k*T[0][1];
    result[0][2] = k*T[0][2];
    result[1][0] = k*T[1][0];
    result[1][1] = k*T[1][1];
    result[1][2] = k*T[1][2];
    result[2][0] = k*T[2][0];
    result[2][1] = k*T[2][1];
    result[2][2] = k*T[2][2];


    return result;
}


// Product of two tensors
Tensor operator*(const Tensor& T1, const Tensor& T2) {

    // Initialize the result
    Tensor R;


    // Compute the product
    R[0][0] = T1[0][0]*T2[0][0] + T1[0][1]*T2[1][0] + T1[0][2]*T2[2][0];
    R[0][1] = T1[0][0]*T2[0][1] + T1[0][1]*T2[1][1] + T1[0][2]*T2[2][1];
    R[0][2] = T1[0][0]*T2[0][2] + T1[0][1]*T2[1][2] + T1[0][2]*T2[2][2];

    R[1][0] = T1[1][0]*T2[0][0] + T1[1][1]*T2[1][0] + T1[1][2]*T2[2][0];
    R[1][1] = T1[1][0]*T2[0][1] + T1[1][1]*T2[1][1] + T1[1][2]*T2[2][1];
    R[1][2] = T1[1][0]*T2[0][2] + T1[1][1]*T2[1][2] + T1[1][2]*T2[2][2];

    R[2][0] = T1[2][0]*T2[0][0] + T1[2][1]*T2[1][0] + T1[2][2]*T2[2][0];
    R[2][1] = T1[2][0]*T2[0][1] + T1[2][1]*T2[1][1] + T1[2][2]*T2[2][1];
    R[2][2] = T1[2][0]*T2[0][2] + T1[2][1]*T2[1][2] + T1[2][2]*T2[2][2];


    return R;
}