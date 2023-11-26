//
// Created by ruben on 26/11/23.
//

#ifndef FLOWBETWEENFLATPLATES_VECTORFIELD_H
#define FLOWBETWEENFLATPLATES_VECTORFIELD_H

#include <vector>
#include "../../math/geometricVector/GeometricVector.h"


class VectorField {
public:
    // Vector field
    std::vector<GeometricVector> field;


    // VectorField constructor and destructor
    VectorField();
    ~VectorField();


    // VectorField methods
    void initialize(int length);
};


#endif //FLOWBETWEENFLATPLATES_VECTORFIELD_H
