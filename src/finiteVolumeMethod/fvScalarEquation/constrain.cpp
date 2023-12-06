//
// Created by ruben on 6/12/23.
//

#include "FvScalarEquation.h"


void FvScalarEquation::constrain(const PolyMesh& theMesh, const ScalarBoundaryConditions& PhiBCs) {

    // Auxiliary variables declaration
    double BCValue;
    int iOwner, iNeighbour, iPeriodicFace, iHalo;
    std::string BCType;


    // Loop over all the boundaries
    for (int k = 0; k < theMesh.nBoundaries; ++k) {

        BCType = PhiBCs.type[k];
        BCValue = PhiBCs.value[k];


        // Loop over all the boundary faces of the k_th boundary
        for (int i = theMesh.boundaries[k].startFace; i < theMesh.boundaries[k].startFace + theMesh.boundaries[k].nBoundaryFaces; ++i) {

            iOwner = theMesh.faces[i].iOwner;
            iNeighbour = theMesh.faces[i].iNeighbour;
            iPeriodicFace = theMesh.faces[i].iPeriodicFace;
            iHalo = theMesh.faces[iPeriodicFace].iOwner;

            if (BCType == "fixedValue") {

                A.addValue(1, {iNeighbour, iNeighbour});
                b.field[iNeighbour] = BCValue;

            } else if (BCType == "zeroGradient" || BCType == "empty") {

                A.addValue(1, {iNeighbour, iNeighbour});
                A.addValue(-1, {iNeighbour, iOwner});
                b.field[iNeighbour] = 0;

            } else if (BCType == "periodic") {

                A.addValue(0.5, {iNeighbour, iNeighbour});
                A.addValue(-0.5, {iNeighbour, iHalo});
                b.field[iNeighbour] = 0;

            } else {

                printf("ERROR. No correct boundary condition type selected !!\n");
            }
        }
    }
}