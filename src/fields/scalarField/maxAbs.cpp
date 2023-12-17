//
// Created by ruben on 7/12/23.
//

#include <cmath>
#include "ScalarField.h"


double ScalarField::maxAbs() {

    double max = 0;

    for (int i = 0; i < size(); ++i) {

        if ( max < fabs(at(i)) ) {

            max = fabs(at(i));
        }
    }


    return max;
}