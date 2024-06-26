//
// Created by ruben on 6/12/23.
//

#include "FvScalarEquation.h"


void FvScalarEquation::constrain(const PolyMesh& theMesh, const ScalarBoundaryConditions& PhiBCs) {

    // Auxiliary variables declaration
    double BCValue, dONMag, SfMag, coeffValue;
    int iOwner, iNeighbour, iPeriodicFace, iHalo;
    std::string BCType;


    // Loop over all the boundaries
    for (int k = 0; k < theMesh.nBoundaries; ++k) {

        BCType = PhiBCs.type[k];
        BCValue = PhiBCs.value[k];


        // Loop over all the boundary faces of the k_th boundary
        for (int i = theMesh.boundaries[k].startFace; i < theMesh.boundaries[k].startFace + theMesh.boundaries[k].nBoundaryFaces; ++i) {

            dONMag = theMesh.faces[i].dONMag;
            SfMag = theMesh.faces[i].SfMag;

            iOwner = theMesh.faces[i].iOwner;
            iPeriodicFace = theMesh.faces[i].iPeriodicFace;
            iHalo = theMesh.faces[iPeriodicFace].iOwner;


            if (BCType == "fixedValue") {

                coeffValue = SfMag/dONMag;

                A.diagValue[iOwner] -= coeffValue;
                b[iOwner] -= BCValue * coeffValue;

            } else if (BCType == "zeroGradient" || BCType == "empty") {

            } else if (BCType == "periodic") {

                coeffValue = SfMag/(2*dONMag);

                A.diagValue[iOwner] -= coeffValue;
                A.addValue(coeffValue, {iOwner, iHalo});

            } else {

                printf("ERROR. No correct boundary condition type selected !!\n");
                std::exit(EXIT_FAILURE);
            }
        }
    }
 }