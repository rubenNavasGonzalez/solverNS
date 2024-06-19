//
// Created by ruben on 7/04/24.
//

#include <iostream>
#include <fstream>
#include "computeEkFull.h"


double computeEkFull(const VectorField& u, const PolyMesh& theMesh, double t) {


    // Initialize Ek variable
    double Ek = 0;


    // Loop over all the interior elements
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        // Update the Ek value
        Ek += theMesh.elements[i].Vf * u[i] * u[i] / 2;
    }


    // Write the kinetic energy to .csv file along with time
    std::ofstream outfileData;
    outfileData.open("Ek.csv", std::ios_base::app);

    if (outfileData.is_open()) {

        // Write the header if file is created
        if (outfileData.tellp() == 0) {

            outfileData << "\"t\"," << "\"Ek\"\n";
        }

        // Append data (time, bulk velocity)
        outfileData << std::to_string(t) << "," << Ek << "\n";
        outfileData.close();
    } else {

        printf("Error. Unable to open file. \n");
        std::exit(EXIT_FAILURE);
    }


    return Ek;
}