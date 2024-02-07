//
// Created by ruben on 31/01/24.
//

#include <fstream>
#include <iostream>
#include "writeTurbulentChannelFlowSetUp.h"


void writeTurbulentChannelFlowSetUp(double delta, double Lx, double Ly, double Lz, double Nx, double Ny, double Nz,
                                    double sx, double sy, double sz, double nu, double yPlusMin, double solverTolerance,
                                    std::string solver, int maxIter, double tInit, double tFinal, double f, double steadyStateCriterion,
                                    double writeIntervalCSV,  double writeIntervalVTK, std::string filename) {


    // Data file
    std::ofstream outFile;
    outFile.open(filename + ".txt");

    if (outFile.is_open()) {

        // Write the header
        outFile << "=== Turbulent Channel Flow Data === \n";
        outFile << "Description: Case set-up data. \n";
        outFile << "=================================== \n\n";


        // Write the geometrical data
        outFile << "Geometrical data: \n";
        outFile << "\tdelta = " << delta << " m\n";
        outFile << "\tLx = " << Lx << " m\n";
        outFile << "\tLy = " << Ly << " m\n";
        outFile << "\tLz = " << Lz << " m\n\n";


        // Write the spatial discretization data
        outFile << "Spatial discretization data: \n";
        outFile << "\tNx = " << Nx << "\n";
        outFile << "\tNy = " << Ny << "\n";
        outFile << "\tNz = " << Nz << "\n";
        outFile << "\tsx = " << sx << "\n";
        outFile << "\tsy = " << sy << "\n";
        outFile << "\tsz = " << sz << "\n";
        outFile << "\tyPlusMin = " << yPlusMin << "\n\n";


        // Write the flow data
        outFile << "Flow data: \n";
        outFile << "\tnu = " << nu << " m^2/s" << "\n\n";


        // Write the linear solver data
        outFile << "Linear solver data: \n";
        outFile << "\tSolver " + solver + "\n";
        outFile << "\tsolverTolerance = " << solverTolerance << "\n";
        outFile << "\tmaxIter = " << maxIter << "\n\n";


        // Write the temporal advancement data
        outFile << "Temporal advancement data: \n";
        outFile << "\ttInit = " << tInit << " s" << "\n";
        outFile << "\ttFinal = " << tFinal << " s" << "\n";
        outFile << "\tf = " << f << "\n";
        outFile << "\tsteadyStateCriterion = " << steadyStateCriterion << "\n\n";


        // Write the file recording data
        outFile << "File recording data: \n";
        outFile << "\twriteInterval2CSV = " << writeIntervalCSV << " s" << "\n";
        outFile << "\twriteInterval2VTK = " << writeIntervalVTK << " s" << "\n\n";

    } else {

        printf("Error. Unable to open file. \n");
    }
}