//
// Created by ruben on 21/01/24.
//

#include "Tensor.h"


// Tensor default constructor and destructor
Tensor::Tensor() {

    this->resize(3, {0,0,0});
}

Tensor::~Tensor() = default;
