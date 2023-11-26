//
// Created by rnava on 11/26/2023.
//

#ifndef FLOWBETWEENFLATPLATES_SCALARBOUNDARYCONDITIONS_H
#define FLOWBETWEENFLATPLATES_SCALARBOUNDARYCONDITIONS_H

#include <vector>
#include <string>


class ScalarBoundaryConditions {
public:
    // Arrays of the boundary conditions of a scalar field (each array position is referred to a boundary)
    std::vector<std::string> type;
    std::vector<double> value;


    // ScalarBoundaryConditions constructor and destructor
    ScalarBoundaryConditions();
    ~ScalarBoundaryConditions();


    // ScalarBoundaryConditions methods
    // Add a new boundary condition
    void addBC(std::string type_, double value_);
};


#endif //FLOWBETWEENFLATPLATES_SCALARBOUNDARYCONDITIONS_H
