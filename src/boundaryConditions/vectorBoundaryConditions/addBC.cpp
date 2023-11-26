//
// Created by ruben on 26/11/23.
//

#include "VectorBoundaryConditions.h"


void VectorBoundaryConditions::addBC(std::string type_, GeometricVector value_) {

    type.push_back(type_);
    value.push_back(value_);
}