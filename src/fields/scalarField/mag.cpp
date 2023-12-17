//
// Created by ruben on 12/12/23.
//

#include <cmath>
#include "ScalarField.h"


double ScalarField::mag() {

    double sum = 0;

    for (int i = 0; i < size(); ++i) {

        sum += pow(at(i),2);
    }


    return sqrt(sum);
}