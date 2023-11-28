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


// Geometric vector operator overloading functions
GeometricVector &GeometricVector::operator=(const GeometricVector &v) {
    this->x = v.x;
    this->y = v.y;
    this->z = v.z;
    return *this;
}

GeometricVector operator+(const GeometricVector &v1, const GeometricVector &v2) {
    return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
}

GeometricVector &GeometricVector::operator+=(const GeometricVector &v) {
    this->x += v.x;
    this->y += v.y;
    this->z += v.z;
    return *this;
}

GeometricVector operator-(const GeometricVector &v1, const GeometricVector &v2) {
    return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
}

GeometricVector &GeometricVector::operator-=(const GeometricVector &v) {
    this->x -= v.x;
    this->y -= v.y;
    this->z -= v.z;
    return *this;
}

double operator*(const GeometricVector &v1, const GeometricVector &v2) {
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

GeometricVector operator*(const double& k, const GeometricVector& v) {
    return { k * v.x, k * v.y, k * v.z };
}

GeometricVector operator*(const GeometricVector &v, const double &k) {
    return { v.x * k, v.y * k, v.z * k };
}

GeometricVector operator/(const GeometricVector &v, const double &k) {
    return { v.x / k, v.y / k, v.z / k };
}

GeometricVector &GeometricVector::operator/=(const double& k) {
    this->x = this->x/k;
    this->y = this->y/k;
    this->z = this->z/k;
}
