//
// Created by rnava on 11/26/2023.
//

#include "ScalarBoundaryConditions.h"


void ScalarBoundaryConditions::addBC(std::string type_, double value_) {

    type.push_back(type_);
    value.push_back(value_);
}