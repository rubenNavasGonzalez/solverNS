//
// Created by ruben on 26/11/23.
//

#ifndef FLOWBETWEENFLATPLATES_VECTORBOUNDARYCONDITIONS_H
#define FLOWBETWEENFLATPLATES_VECTORBOUNDARYCONDITIONS_H

#include <vector>
#include <string>
#include "../../math/geometricVector/GeometricVector.h"


class VectorBoundaryConditions {
public:
    // Arrays of the boundary conditions of a vector field (each array position is referred to a boundary)
    std::vector<std::string> type;
    std::vector<GeometricVector> value;


    // VectorBoundaryConditions constructor and destructor
    VectorBoundaryConditions();
    ~VectorBoundaryConditions();


    // VectorBoundaryConditions methods
    // Add a new boundary condition
    void addBC(std::string type_, GeometricVector value_);
};


#endif //FLOWBETWEENFLATPLATES_VECTORBOUNDARYCONDITIONS_H
