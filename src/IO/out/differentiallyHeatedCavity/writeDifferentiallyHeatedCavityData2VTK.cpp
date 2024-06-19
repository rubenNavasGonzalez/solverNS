//
// Created by ruben on 6/04/24.
//

#include <iostream>
#include <fstream>
#include "writeDifferentiallyHeatedCavityData2VTK.h"

void writeDifferentiallyHeatedCavityData2VTK(const PolyMesh& theMesh, const ScalarField& p, const ScalarField& T, const VectorField& u,
                                       const VectorField& omega, const ScalarField& nut, const ScalarField& QCrit,
                                       const ScalarBoundaryConditions& pBCs, const ScalarBoundaryConditions& TBCs,
                                       const VectorBoundaryConditions& uBCs, const VectorBoundaryConditions& omegaBCs,
                                       const ScalarBoundaryConditions& nutBCs, const ScalarBoundaryConditions& QCritBCs, double t) {


    // Auxiliary variables declaration
    std::string filename = "Time_" + std::to_string(t);


    // Write the mesh data
    theMesh.writeMesh2VTK(filename);


    // Specify the number of data units
    std::ofstream outfile;
    outfile.open(filename + ".vtk", std::ios_base::app);

    if (outfile.is_open()) {

        outfile << "POINT_DATA " << theMesh.nNodes << "\n\n";
        outfile.close();
    } else {

        printf("Unable to open file\n");
        std::exit(EXIT_FAILURE);
    }


    // Write the pressure data
    p.writeScalarField2VTK(filename, "p", theMesh, pBCs);

    // Write the temperature data
    T.writeScalarField2VTK(filename, "T", theMesh, TBCs);


    // Write the QCriterion data
    QCrit.writeScalarField2VTK(filename, "QCrit", theMesh, QCritBCs);


    // Write the turbulent viscosity data
    nut.writeScalarField2VTK(filename, "nut", theMesh, nutBCs);


    // Write the velocity data
    u.writeVectorField2VTK(filename, "U", theMesh, uBCs);


    // Write the vorticity data
    omega.writeVectorField2VTK(filename, "omega", theMesh, omegaBCs);
}