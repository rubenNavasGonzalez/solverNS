//
// Created by ruben on 7/12/23.
//

#include <cmath>
#include "ScalarField.h"


double ScalarField::maxAbs() {

    double max = 0;

    for (int i = 0; i < field.size(); ++i) {

        if ( max < fabs(field[i]) ) {

            max = fabs(field[i]);
        }
    }

    return max;
}