//
// Created by ruben on 6/02/24.
//

#include <string>
#include "initializeFromSelectedTime.h"
#include "../../IO/in/readCSV.h"


VectorField initializeFromSelectedTime(const PolyMesh& theMesh, double t) {

    // Initialize the velocity field
    VectorField u;
    u.assign(theMesh.nInteriorElements, {0,0,0});


    // Read the .csv file according to the selected time
    std::vector < std::vector<double> > dataCSV = readCSV("Time_" + std::to_string(t) + ".csv");


    // Loop over all the interior elements and assign the read velocity to the new velocity field
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        u[i].x = dataCSV[i][4];
        u[i].y = dataCSV[i][5];
        u[i].z = dataCSV[i][6];
    }

    return u;
}