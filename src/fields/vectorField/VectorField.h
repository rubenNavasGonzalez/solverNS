//
// Created by ruben on 26/11/23.
//

#ifndef FLOWBETWEENFLATPLATES_VECTORFIELD_H
#define FLOWBETWEENFLATPLATES_VECTORFIELD_H

#include <vector>
#include "../../math/geometricVector/GeometricVector.h"
#include "../../mesh/polyMesh/PolyMesh.h"
#include "../../boundaryConditions/vectorBoundaryConditions/VectorBoundaryConditions.h"


class VectorField : public std::vector<GeometricVector> {
public:

    // VectorField constructor and destructor
    VectorField();
    ~VectorField();


    // VectorField methods
    double maxAbs();
    void writeVectorField2VTK(const std::string& filename, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs);
    friend VectorField operator+(const VectorField& field1, const VectorField& field2);
    friend VectorField operator-(const VectorField& field1, const VectorField& field2);
    friend VectorField operator*(const double& k, const VectorField& field);
    friend VectorField operator/(const VectorField& field, const double& k);
};


#endif //FLOWBETWEENFLATPLATES_VECTORFIELD_H
