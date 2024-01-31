//
// Created by ruben on 31/01/24.
//

#include "Tensor.h"


Tensor Tensor::skewSymmetric() {

    // Initialize the skew-symmetric part
    Tensor S;


    // Compute the skew-symmetric part
    S = 0.5*(*this - this->transpose());


    return S;
}