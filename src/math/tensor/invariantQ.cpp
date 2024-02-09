//
// Created by ruben on 8/02/24.
//

#include <cmath>
#include "Tensor.h"


double Tensor::invariantQ() {

    return 0.5*(pow(this->trace(),2) - ((*this) * (*this)).trace());
}