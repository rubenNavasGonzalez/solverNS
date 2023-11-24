//
// Created by ruben on 19/11/23.
//

#include <cmath>
#include "GeometricVector.h"


// GeometricVector default constructor and destructor
GeometricVector::GeometricVector() = default;
GeometricVector::~GeometricVector() = default;

// GeometricVector constructor from a coordinate
GeometricVector::GeometricVector(double x_, double y_, double z_) {
    x = x_;
    y = y_;
    z = z_;
}


// GeometricVector magnitude
double GeometricVector::mag() {
    return sqrt(x*x + y*y + z*z);
}
