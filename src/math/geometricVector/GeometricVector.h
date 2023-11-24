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
};


#endif //FLOWBETWEENFLATPLATES_GEOMETRICVECTOR_H
