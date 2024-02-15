//
// Created by ruben on 8/02/24.
//

#include "Tensor.h"


double Tensor::invariantV2() {

    // Auxiliary variables definition
    Tensor S = this->symmetric();
    Tensor W = this->skewSymmetric();


    return 4*(( (S*S) * (W*W) ).trace() - 2*S.invariantQ()*W.invariantQ());
}