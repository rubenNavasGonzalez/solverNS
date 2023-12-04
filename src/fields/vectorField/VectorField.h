//
// Created by ruben on 26/11/23.
//

#ifndef FLOWBETWEENFLATPLATES_VECTORFIELD_H
#define FLOWBETWEENFLATPLATES_VECTORFIELD_H

#include <vector>
#include "../../math/geometricVector/GeometricVector.h"
#include "../../mesh/polyMesh/PolyMesh.h"
#include "../../boundaryConditions/vectorBoundaryConditions/VectorBoundaryConditions.h"


class VectorField {
public:
    // Vector field
    std::vector<GeometricVector> field;


    // VectorField constructor and destructor
    VectorField();
    ~VectorField();


    // VectorField methods
    void initialize(int length);
    void applyBCs(const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs);
    friend VectorField operator+(const VectorField& field1, const VectorField& field2);
    friend VectorField operator-(const VectorField& field1, const VectorField& field2);
};


#endif //FLOWBETWEENFLATPLATES_VECTORFIELD_H
