//
// Created by ruben on 8/02/24.
//

#include "Tensor.h"


double Tensor::trace() {

    return this->at(0).at(0) + this->at(1).at(1) + this->at(2).at(2);
}