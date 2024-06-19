//
// Created by ruben on 7/04/24.
//

#include <iostream>
#include <fstream>
#include "writeDifferentiallyHeatedCavityResults2CSV.h"


void writeDifferentiallyHeatedCavityResults2CSV(const PolyMesh& theMesh, const VectorField& u, const ScalarField& p, const ScalarField& T,
                                                const ScalarField& nut, double t) {

    // Data file
    std::ofstream outfileData;
    outfileData.open("Time_" + std::to_string(t) + ".csv");

    if (outfileData.is_open()) {

        // Write the data file header
        outfileData << "\"x\"," << "\"y\"," << "\"z\"," << "\"p\"," << "\"uX\"," << "\"uY\"," << "\"uZ\"," << "\"T\","
                    << "\"nut\"\n";

        // Assemble the .csv file
        for (int i = 0; i < p.size(); ++i) {

            outfileData << theMesh.elements[i].centroid.x << "," << theMesh.elements[i].centroid.y << "," << theMesh.elements[i].centroid.z << ","
                        << p[i] << "," << u[i].x << "," << u[i].y << "," << u[i].z << "," << T[i] << "," << nut[i] << "\n";
        }

        // Close the file
        outfileData.close();

    } else {

        printf("Error. Unable to open file. \n");
        std::exit(EXIT_FAILURE);
    }


    // File with the write time
    std::ofstream outfileTime;
    outfileTime.open("writeTime.txt", std::ios_base::app);

    if (outfileTime.is_open()) {

        outfileTime << std::to_string(t) << "\n";
        outfileTime.close();
    } else {

        printf("Error. Unable to open file. \n");
        std::exit(EXIT_FAILURE);
    }
}