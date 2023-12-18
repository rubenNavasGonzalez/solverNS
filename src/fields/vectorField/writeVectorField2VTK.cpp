//
// Created by ruben on 18/12/23.
//

#include <fstream>
#include <iostream>
#include "VectorField.h"
#include "../../interpolation/interpolateFromElements2Nodes/interpolateVectorFieldFromElements2Nodes.h"


void VectorField::writeVectorField2VTK(const std::string& filename, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs) {

    // Write the Mesh into .VTK format
    theMesh.writeMesh2VTK(filename);


    // Open the mesh file in append mode
    std::ofstream outfile;
    outfile.open(filename + ".vtk", std::ios_base::app);


    // Interpolate the values from the elements to the nodes
    VectorField nodeField = interpolateVectorFieldFromElements2Nodes(*this, theMesh, PhiBCs);


    // Write cell based and point based fields
    if (outfile.is_open()) {
        // Write the scalar field (cell based)
        /*outfile << "CELL_DATA " << mesh.nElements << "\n";
        outfile << "VECTORS vectorField float" << "\n";
        for (int i = 0; i < mesh.nElements; i++) {
            outfile << field.field[i] << "\n";
        }*/

        // Write the scalar field (point based)
        outfile << "POINT_DATA " << theMesh.nNodes << "\n";
        outfile << "VECTORS vectorField float" << "\n";
        for (int i = 0; i < theMesh.nNodes; i++) {
            outfile << nodeField[i].x << "\t" << nodeField[i].y << "\t" << nodeField[i].z << "\n";
        }

        // Close the file
        outfile.close();

    } else {
        std::cout << "Unable to open file";
    }
}