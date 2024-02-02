//
// Created by ruben on 21/01/24.
//

#include <iostream>
#include <fstream>
#include "writeTurbulentChannelFlowData2CSV.h"


void writeTurbulentChannelFlowData2CSV(const PolyMesh& theMesh, const ScalarField& p, const VectorField& u, const VectorField& omega,
                                       const ScalarField& nut, const TensorField& gradU, double uBulk, double t) {

    // Data file
    std::ofstream outfileData;
    outfileData.open("Time_" + std::to_string(t) + ".csv");

    if (outfileData.is_open()) {


        // Write the data file header
        outfileData << "x," << "y," << "z," << "p," << "uX," << "uY," << "uZ," << "wX," << "wY," << "wZ," << "nut," <<
        "partU_partX," << "partV_partX," << "partW_partX," << "partU_partY," << "partV_partY," << "partW_partY," <<
        "partU_partZ," << "partV_partZ," << "partW_partZ," << "uBulk \n";


        // Assemble the .csv file
        for (int i = 0; i < p.size(); ++i) {

            outfileData << theMesh.elements[i].centroid.x << "," << theMesh.elements[i].centroid.y << "," << theMesh.elements[i].centroid.z << "," <<
            p[i] << "," << u[i].x << "," << u[i].y << "," << u[i].z << "," << omega[i].x << "," << omega[i].y << "," << omega[i].z << "," <<
            nut[i] << "," << gradU[i][0][0]<< "," << gradU[i][0][1]<< "," << gradU[i][0][2]<< "," << gradU[i][1][0]<< "," << gradU[i][1][1] << "," <<
            "," << gradU[i][1][2] <<"," << gradU[i][2][0] << "," << gradU[i][2][1] << "," << gradU[i][2][2] << "," << uBulk << "\n";
        }


        // Close the file
        outfileData.close();

    } else {

        printf("Error. Unable to open file. \n");
    }


    // File with the write time
    std::ofstream outfileTime;
    outfileTime.open("writeTime.txt", std::ios_base::app);

    if (outfileTime.is_open()) {

        outfileTime << std::to_string(t) << "\n";
        outfileTime.close();
    } else {

        printf("Error. Unable to open file. \n");
    }
}