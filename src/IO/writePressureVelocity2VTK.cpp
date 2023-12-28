//
// Created by ruben on 28/12/23.
//

#include <fstream>
#include <iostream>
#include <string>
#include "writePressureVelocity2VTK.h"
#include "../interpolation/interpolateFromElements2Nodes/interpolateScalarFieldFromElements2Nodes.h"
#include "../interpolation/interpolateFromElements2Nodes/interpolateVectorFieldFromElements2Nodes.h"


void writePressureVelocity2VTK(const PolyMesh& theMesh, const ScalarField& p, const VectorField& u, const ScalarBoundaryConditions& pBCs,
                               const VectorBoundaryConditions& uBCs, double t) {

    std::ofstream outfile;
    outfile.open("Time_" + std::to_string(t) + ".vtk");

    if (outfile.is_open()) {
        // Write the header
        outfile << "# vtk DataFile Version 3.0" << "\n";
        outfile << "Unstructured Mesh" << "\n";
        outfile << "ASCII" << "\n\n";

        // Write the points
        int nPoints = theMesh.nodes.size();
        outfile << "DATASET UNSTRUCTURED_GRID" << "\n";
        outfile << "POINTS " << nPoints << " double" << "\n\n";
        for (int i = 0; i < nPoints; i++) {
            outfile << theMesh.nodes[i].x << " " << theMesh.nodes[i].y << " " << theMesh.nodes[i].z << "\n";
        }

        // Write the cells
        int nCells = theMesh.elements.size();
        int cellType = 8;
        int num_total_cells = nCells * (1 + cellType);
        outfile << "CELLS " << nCells << " " << num_total_cells << "\n";
        for (int i = 0; i < nCells; i++) {
            outfile << cellType;
            for (int j = 0; j < cellType; j++) {
                outfile << " " << theMesh.elements[i].iNodes[j];
            }
            outfile << std::endl;
        }

        // Write the cell types
        outfile << "CELL_TYPES " << nCells << "\n";
        for (int i = 0; i < nCells; i++) {
            outfile << theMesh.elements[i].elementType << "\n";
        }
        outfile << "\n";

        // Close the file
        outfile.close();
    } else {
        std::cout << "Unable to open file";
    }


    outfile.open("Time_" + std::to_string(t) + ".vtk", std::ios_base::app);

    // Interpolate the values from the elements to the nodes
    ScalarField nodeField = interpolateScalarFieldFromElements2Nodes(p, theMesh, pBCs);


    // Write cell based and point based fields
    if (outfile.is_open()) {
        // Write the scalar field (cell based)
        /*outfile << "CELL_DATA " << mesh.nElements << "\n";
        outfile << "SCALARS " + filename + " float 1" << "\n";
        outfile << "LOOKUP_TABLE default" << "\n";
        for (int i = 0; i < mesh.nElements; i++) {
            outfile << field.field[i] << "\n";
        }*/

        // Write the scalar field (point based)
        outfile << "POINT_DATA " << theMesh.nNodes << "\n\n";
        outfile << "SCALARS p float 1" << "\n";
        outfile << "LOOKUP_TABLE default" << "\n";
        for (int i = 0; i < theMesh.nNodes; i++) {
            outfile << nodeField[i] << "\n";
        }
        outfile << "\n";

        // Close the file
        outfile.close();

    } else {
        std::cout << "Unable to open file";
    }


    outfile.open("Time_" + std::to_string(t) + ".vtk", std::ios_base::app);

    // Interpolate the values from the elements to the nodes
    VectorField uNodeField = interpolateVectorFieldFromElements2Nodes(u, theMesh, uBCs);


    // Write cell based and point based fields
    if (outfile.is_open()) {
        // Write the scalar field (cell based)
        /*outfile << "CELL_DATA " << mesh.nElements << "\n";
        outfile << "VECTORS " + filename + " float" << "\n";
        for (int i = 0; i < mesh.nElements; i++) {
            outfile << field.field[i] << "\n";
        }*/

        // Write the scalar field (point based)
        outfile << "VECTORS U float" << "\n";
        for (int i = 0; i < theMesh.nNodes; i++) {
            outfile << uNodeField[i].x << "\t" << uNodeField[i].y << "\t" << uNodeField[i].z << "\n";
        }

        // Close the file
        outfile.close();

    } else {
        std::cout << "Unable to open file";
    }
}