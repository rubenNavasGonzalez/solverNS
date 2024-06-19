//
// Created by ruben on 10/04/24.
//

#include "initializeFromSelectedTime.h"
#include "../../IO/in/readCSV.h"


ScalarField initializeScalarFieldFromSelectedTime(const PolyMesh& theMesh, double t) {

    // Initialize the velocity field
    ScalarField T;
    T.assign(theMesh.nInteriorElements, 0);


    // Read the .csv file according to the selected time
    std::vector < std::vector<double> > dataCSV = readCSV("Time_" + std::to_string(t) + ".csv");


    // Loop over all the interior elements and assign the read velocity to the new velocity field
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        T[i] = dataCSV[i][7];
    }

    return T;
}