//
// Created by ruben on 20/12/23.
//

#include <cmath>
#include "VectorField.h"
#include "../scalarField/ScalarField.h"


double VectorField::maxAbs() {

    double xAbs, yAbs, zAbs;

    ScalarField maxAbsComponent;
    maxAbsComponent.assign(size(), 0);

    for (int i = 0; i < size(); ++i) {

        xAbs = this->at(i).x;
        yAbs = this->at(i).y;
        zAbs = this->at(i).z;

        maxAbsComponent[i] = fmax(fmax(xAbs, yAbs), zAbs);
    }


    return maxAbsComponent.max();
}