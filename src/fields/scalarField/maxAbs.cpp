//
// Created by ruben on 7/12/23.
//

#include <cmath>
#include "ScalarField.h"


double ScalarField::maxAbs() {

    // Initialize the result
    double max = 0;


    // Perform the operation
    for (int i = 0; i < this->size(); ++i) {

        if ( max < fabs(this->at(i)) ) {

            max = fabs(this->at(i));
        }
    }


    return max;
}