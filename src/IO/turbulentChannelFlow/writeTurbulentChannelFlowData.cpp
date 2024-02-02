//
// Created by ruben on 31/01/24.
//

#include <fstream>
#include <iostream>
#include "writeTurbulentChannelFlowData.h"


void writeTurbulentChannelFlowData(double delta, double Lx, double Ly, double Lz, double Nx, double Ny, double Nz,
                                   double sx, double sy, double sz, double nu, double yPlusMin, double solverTolerance,
                                   std::string solver, int maxIter, double tInit, double tFinal, double f,
                                   double steadyStateCriterion, int writeInterval) {

    // Data file
    std::ofstream outFile;
    outFile.open("caseSetUp.txt");

    if (outFile.is_open()) {

        // Write the header
        outFile << "=== Turbulent Channel Flow Data === \n";
        outFile << "Description: Case set-up data. \n";
        outFile << "=================================== \n\n";


        // Write the geometrical data
        outFile << "Geometrical data: \n";
        outFile << "\tdelta = " << delta << "\n";
        outFile << "\tLx = " << Lx << "\n";
        outFile << "\tLy = " << Ly << "\n";
        outFile << "\tLz = " << Lz << "\n\n";


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
        outFile << "\tnu = " << nu << "\n\n";


        // Write the linear solver data
        outFile << "Linear solver data: \n";
        outFile << "\tSolver " + solver + ".\n";
        outFile << "\tsolverTolerance = " << solverTolerance << "\n";
        outFile << "\tmaxIter = " << maxIter << "\n\n";


        // Write the temporal advancement data
        outFile << "Temporal advancement data: \n";
        outFile << "\ttInit = " << tInit << "\n";
        outFile << "\ttFinal = " << tFinal << "\n";
        outFile << "\tf = " << f << "\n";
        outFile << "\tsteadyStateCriterion = " << steadyStateCriterion << "\n\n";


        // Write the file recording data
        outFile << "File recording data: \n";
        outFile << "\twriteInterval = " << writeInterval << "\n\n";

    } else {

        printf("Error. Unable to open file. \n");
    }
}