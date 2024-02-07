//
// Created by ruben on 7/02/24.
//

#ifndef SOLVERNS_PROBE_H
#define SOLVERNS_PROBE_H

#include <vector>
#include <string>
#include "../../fields/vectorField/VectorField.h"


class Probe {
public:
    // Probe attributes
    Node p;                                 // Probe coordinates
    std::vector<int> iElements;             // Assigned elements (closest two)
    std::vector<double> d;                  // Distance of the assigned elements to the probe


    // Probe constructor and destructor
    Probe();
    ~Probe();
    Probe(double x, double y, double z);


    // Probe methods
    void assign(const PolyMesh& theMesh);
    void writeField(const VectorField& Phi, double t, std::string filename);
};


#endif //SOLVERNS_PROBE_H
