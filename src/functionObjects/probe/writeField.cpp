//
// Created by ruben on 7/02/24.
//

#include <fstream>
#include "Probe.h"


void Probe::writeField(const VectorField& Phi, double t, std::string filename) {

    // Auxiliary variable definition
    double dTotal;
    GeometricVector PhiProbe;


    // Compute the Phi variable at probe position
    dTotal = this->d[0] + this->d[1];
    PhiProbe = Phi[this->iElements[0]]*(1 - this->d[0]/dTotal) + Phi[this->iElements[1]]*(1 - this->d[1]/dTotal);


    // Write the probe data to .csv file along with time
    std::ofstream outfileData;
    outfileData.open(filename + ".csv", std::ios_base::app);

    if (outfileData.is_open()) {

        // Write the header if file is created
        if (outfileData.tellp() == 0) {

            outfileData << "\"t\"," << "\"probeX\"," << "\"probeY\"," << "\"probeZ\"\n";
        }

        // Append data (time, bulk velocity)
        outfileData << std::to_string(t) << "," << PhiProbe.x << "," << PhiProbe.y << "," << PhiProbe.z << "\n";
        outfileData.close();
    } else {

        printf("Error. Unable to open file. \n");
        std::exit(EXIT_FAILURE);
    }
}