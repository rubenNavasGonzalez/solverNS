//
// Created by ruben on 31/01/24.
//

#include "Tensor.h"


Tensor Tensor::symmetric() {

    // Initialize the symmetric part
    Tensor S;


    // Compute the symmetric part
    S = 0.5*(*this + this->transpose());


    return S;
}