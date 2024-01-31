//
// Created by ruben on 26/11/23.
//

#ifndef FLOWBETWEENFLATPLATES_SCALARFIELD_H
#define FLOWBETWEENFLATPLATES_SCALARFIELD_H

#include <vector>
#include <string>
#include "../../boundaryConditions/scalarBoundaryConditions/ScalarBoundaryConditions.h"
#include "../../mesh/polyMesh/PolyMesh.h"


class ScalarField : public std::vector<double> {
public:
    // ScalarField constructor and destructor
    ScalarField();
    ~ScalarField();


    // ScalarField methods
    double max();
    double min();
    double maxAbs();
    double rms();
    void writeScalarField2VTK(const std::string& filename, const std::string& field, const PolyMesh& theMesh, const ScalarBoundaryConditions& PhiBCs) const;
    friend ScalarField operator-(const ScalarField& field1, const ScalarField& field2);
    friend ScalarField operator+(const ScalarField& field1, const ScalarField& field2);
    friend double operator*(const ScalarField& field1, const ScalarField& field2);
    friend ScalarField operator*(const double& k, const ScalarField& field);
    friend ScalarField operator/(const ScalarField& field, const double& k);
};


#endif //FLOWBETWEENFLATPLATES_SCALARFIELD_H
