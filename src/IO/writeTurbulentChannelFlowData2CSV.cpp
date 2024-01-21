//
// Created by ruben on 21/01/24.
//

#include <iostream>
#include <fstream>
#include "writeTurbulentChannelFlowData2CSV.h"


void writeTurbulentChannelFlowData2CSV(const ScalarField& p, const VectorField& u, const VectorField& omega, double uBulk, double t) {

    std::ofstream outfile;
    outfile.open("Time_" + std::to_string(t) + ".txt");

    if (outfile.is_open()) {

        // Write the data file header
        outfile << "p," << "uX," << "uY," << "uZ," << "wX," << "wY," << "wZ," << "uBulk \n";

        // Assemble the .csv file
        for (int i = 0; i < p.size(); ++i) {

            outfile << p[i] << "," << u[i].x << "," << u[i].y << "," << u[i].z << "," << omega[i].x << "," << omega[i].y << "," << omega[i].z << "," << uBulk << "\n";
        }

    } else {

        printf("Error. Unable to open file. \n");
    }
}