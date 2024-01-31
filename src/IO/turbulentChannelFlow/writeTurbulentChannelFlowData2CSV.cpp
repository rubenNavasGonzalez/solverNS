//
// Created by ruben on 21/01/24.
//

#include <iostream>
#include <fstream>
#include "writeTurbulentChannelFlowData2CSV.h"


void writeTurbulentChannelFlowData2CSV(const PolyMesh& theMesh, const ScalarField& p, const VectorField& u, const VectorField& omega,
                                       const ScalarField& nut, double uBulk, double t) {

    // Data file
    std::ofstream outfileData;
    outfileData.open("Time_" + std::to_string(t) + ".csv");

    if (outfileData.is_open()) {

        // Write the data file header
        outfileData << "x," << "y," << "z," << "p," << "uX," << "uY," << "uZ," << "wX," << "wY," << "wZ," << "nut," << "uBulk \n";

        // Assemble the .csv file
        for (int i = 0; i < p.size(); ++i) {

            outfileData << theMesh.elements[i].centroid.x << "," << theMesh.elements[i].centroid.y << "," << theMesh.elements[i].centroid.z << "," << p[i] << "," << u[i].x << "," << u[i].y << "," << u[i].z << "," << omega[i].x << "," << omega[i].y << "," << omega[i].z << "," << nut[i] << "," << uBulk << "\n";
        }

    } else {

        printf("Error. Unable to open file. \n");
    }


    // File with the write time
    std::ofstream outfileTime;
    outfileTime.open("writeTime.txt", std::ios_base::app);

    if (outfileTime.is_open()) {

        outfileTime << std::to_string(t) << "\n";
    } else {

        printf("Error. Unable to open file. \n");
    }
}