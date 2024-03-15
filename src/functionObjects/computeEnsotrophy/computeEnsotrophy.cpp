//
// Created by ruben on 6/03/24.
//

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "computeEnsotrophy.h"


double computeEnsotrophy(const VectorField& omega, const PolyMesh& theMesh, double t) {


    // Initialize Dzeta variable
    double Dzeta = 0;


    // Loop over all the interior elements
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        // Update the Ek value
        Dzeta += omega[i] * omega[i];
    }


    // Average the Ek value by the volume
    Dzeta /= (2*theMesh.nInteriorElements);


    // Write the kinetic energy to .csv file along with time
    std::ofstream outfileData;
    outfileData.open("Ensotrophy.csv", std::ios_base::app);

    if (outfileData.is_open()) {

        // Write the header if file is created
        if (outfileData.tellp() == 0) {

            outfileData << "\"t\"," << "\"Dzeta\"\n";
        }

        // Append data (time, bulk velocity)
        outfileData << std::to_string(t) << "," << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1) << Dzeta << "\n";
        outfileData.close();
    } else {

        printf("Error. Unable to open file. \n");
        std::exit(EXIT_FAILURE);
    }


    return Dzeta;
}