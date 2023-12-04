//
// Created by ruben on 19/11/23.
//

#ifndef FLOWBETWEENFLATPLATES_GEOMETRICVECTOR_H
#define FLOWBETWEENFLATPLATES_GEOMETRICVECTOR_H


class GeometricVector {
public:
    // Vector x,y,z coordinates
    double x{}, y{}, z{};


    // GeometricVector constructors and destructor
    GeometricVector();
    ~GeometricVector();
    GeometricVector(double x_, double y_, double z_);


    // GeometricVector methods
    double mag();
    GeometricVector& operator=(const GeometricVector& v);
    friend GeometricVector operator+(const GeometricVector& v1, const GeometricVector& v2);
    GeometricVector& operator+=(const GeometricVector& v);
    friend GeometricVector operator-(const GeometricVector& v1, const GeometricVector& v2);
    GeometricVector& operator-=(const GeometricVector& v);
    friend double operator*(const GeometricVector& v1, const GeometricVector& v2);
    friend GeometricVector operator*(const double& k, const GeometricVector& v);
    friend GeometricVector operator*(const GeometricVector& v, const double& k);
    friend GeometricVector operator/(const GeometricVector& v, const double& k);
    friend GeometricVector operator/(const GeometricVector& v1, const GeometricVector& v2);
    GeometricVector& operator/=(const double& k);
};


#endif //FLOWBETWEENFLATPLATES_GEOMETRICVECTOR_H
