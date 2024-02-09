//
// Created by ruben on 8/02/24.
//

#include "Tensor.h"


double Tensor::invariantV2() {

    // Auxiliary variables definition
    Tensor S = this->symmetric();
    Tensor Omega = this->skewSymmetric();


    return 4*(( (S*S) * (Omega*Omega) ).trace() - 2*S.invariantQ()*Omega.invariantQ());
}