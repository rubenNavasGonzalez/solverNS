//
// Created by ruben on 31/01/24.
//

#include "Tensor.h"


Tensor Tensor::transpose() {

    // Tensor initialization
    Tensor T;


    // Assembly the transposed tensor
    T[0][0] = this->at(0).at(0);
    T[0][1] = this->at(1).at(0);
    T[0][2] = this->at(2).at(0);
    T[1][0] = this->at(0).at(1);
    T[1][1] = this->at(1).at(1);
    T[1][2] = this->at(2).at(1);
    T[2][0] = this->at(0).at(2);
    T[2][1] = this->at(1).at(2);
    T[2][2] = this->at(2).at(2);


    return T;
}