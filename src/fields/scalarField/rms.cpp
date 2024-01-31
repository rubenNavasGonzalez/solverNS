//
// Created by ruben on 31/01/24.
//

#include <cmath>
#include "ScalarField.h"


double ScalarField::rms() {

    // rms initialization
    double rms = 0;


    // Loop for all the vector elements and compute the quadratic contribution
    for (int i = 0; i < this->size(); ++i) {

        rms += pow(this->at(i), 2);
    }


    return sqrt(rms/this->size());
}