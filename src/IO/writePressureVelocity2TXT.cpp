//
// Created by ruben on 28/12/23.
//

#include <fstream>
#include <iostream>
#include "writePressureVelocity2TXT.h"


void writePressureVelocity2TXT(const PolyMesh& theMesh, const ScalarField& p, const VectorField& u, double t) {

    std::ofstream outfile;
    outfile.open("Time_" + std::to_string(t) + ".txt");

    if (outfile.is_open()) {
        // Write the header
        outfile << "Velocity_x - Velocity_y - Velocity_z - Pressure\n";

        // Write the data of each interior element
        for (int i = 0; i < theMesh.nInteriorElements; ++i) {

            outfile << u[i].x << "\t" << u[i].y << "\t" << u[i].z << "\t" << p[i] << "\n";
        }

        // Close the file
        outfile.close();
    } else {
        std::cout << "Unable to open file";
    }
}