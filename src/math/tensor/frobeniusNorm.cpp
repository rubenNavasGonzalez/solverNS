//
// Created by ruben on 31/01/24.
//

#include <cmath>
#include "Tensor.h"


double Tensor::frobeniusNorm() {

    // Frobenius norm initialization
    double F = 0;


    // Compute the Frobenius norm
    // Loop for all the tensor terms
    for (int i = 0; i < this->size(); ++i) {

        for (int j = 0; j < this[0].size(); ++j) {

            F += pow(this->at(i).at(j), 2);
        }
    }


    return sqrt(F);
}