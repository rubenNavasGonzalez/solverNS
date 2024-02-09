//
// Created by ruben on 25/11/23.
//

#include <fstream>
#include <iostream>
#include "PolyMesh.h"


void PolyMesh::writeMesh2VTK(const std::string& filename) const {

    std::ofstream outfile;
    outfile.open(filename + ".vtk");

    if (outfile.is_open()) {
        // Write the header
        outfile << "# vtk DataFile Version 3.0" << "\n";
        outfile << "Unstructured Mesh" << "\n";
        outfile << "ASCII" << "\n";

        // Write the points
        int nPoints = nodes.size();
        outfile << "DATASET UNSTRUCTURED_GRID" << "\n";
        outfile << "POINTS " << nPoints << " double" << "\n";
        for (int i = 0; i < nPoints; i++) {
            outfile << nodes[i].x << " " << nodes[i].y << " " << nodes[i].z << "\n";
        }

        // Write the cells
        int nCells = elements.size();
        int cellType = 8;
        int num_total_cells = nCells * (1 + cellType);
        outfile << "CELLS " << nCells << " " << num_total_cells << "\n";
        for (int i = 0; i < nCells; i++) {
            outfile << cellType;
            for (int j = 0; j < cellType; j++) {
                outfile << " " << elements[i].iNodes[j];
            }
            outfile << std::endl;
        }

        // Write the cell types
        outfile << "CELL_TYPES " << nCells << "\n";
        for (int i = 0; i < nCells; i++) {
            outfile << elements[i].elementType << "\n";
        }
        outfile << "\n";

        // Close the file
        outfile.close();
    } else {

        std::cout << "Unable to open file\n";
        std::exit(EXIT_FAILURE);
    }
}