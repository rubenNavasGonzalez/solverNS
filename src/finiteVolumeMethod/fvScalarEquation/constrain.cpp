//
// Created by ruben on 6/12/23.
//

#include "FvScalarEquation.h"


void FvScalarEquation::constrain(const PolyMesh& theMesh, double kValue, const ScalarBoundaryConditions& PhiBCs) {

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
            iNeighbour = theMesh.faces[i].iNeighbour;
            iPeriodicFace = theMesh.faces[i].iPeriodicFace;
            iHalo = theMesh.faces[iPeriodicFace].iOwner;


            if (BCType == "fixedValue") {

                coeffValue = SfMag/dONMag;

                A.diagValue[iOwner] -= coeffValue;
                b[iNeighbour] -= BCValue*kValue;

            } else if (BCType == "zeroGradient" || BCType == "empty") {

            } else if (BCType == "periodic") {

                coeffValue = SfMag/(2*dONMag);

                A.diagValue[iOwner] -= coeffValue;
                A.addValue(coeffValue, {iOwner, iHalo});

            } else {

                printf("ERROR. No correct boundary condition type selected !!\n");
            }
        }
    }

    A.diagValue[0] = A.diagValue[0]*1.1;

    for (int i = 0; i < A.diagValue.size(); ++i) {
        A.diagValue[i] = A.diagValue[i]*(-1);
    }

    for (int i = 0; i < A.lowerValue.size(); ++i) {
        A.lowerValue[i] = A.lowerValue[i]*(-1);
    }

    for (int i = 0; i < A.lowerValue.size(); ++i) {
        A.upperValue[i] = A.upperValue[i]*(-1);
    }

    for (int i = 0; i < b.size(); ++i) {
        b[i] = b[i]*(-1);
    }
 }