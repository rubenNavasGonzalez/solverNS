//
// Created by ruben on 3/12/23.
//

#include "VectorField.h"


/*
void VectorField::applyBCs(const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs) {

    // Variables preallocation
    int iOwner, iNeighbour, iHalo, iPeriodicFace;
    GeometricVector PhiOwner, PhiNeighbour, PhiHalo, DivPhiF, BCValue;
    std::string BCType;


    // Loop over all the boundaries
    for (int k = 0; k < theMesh.nBoundaries; ++k) {

        BCType = PhiBCs.type[k];
        BCValue = PhiBCs.value[k];


        // Loop over all the boundary faces of the k_th boundary
        for (int i = theMesh.boundaries[k].startFace; i < theMesh.boundaries[k].startFace + theMesh.boundaries[k].nBoundaryFaces; ++i) {

            iOwner = theMesh.faces[i].iOwner;
            iNeighbour = theMesh.faces[i].iNeighbour;

            if (BCType == "fixedValue") {

                this->field[iNeighbour] = BCValue;
            } else if (BCType == "zeroGradient") {

                this->field[iNeighbour] = this->field[iOwner];
            } else if (BCType == "periodic") {

                iPeriodicFace = theMesh.faces[i].iPeriodicFace;
                iHalo = theMesh.faces[iPeriodicFace].iOwner;

                PhiOwner = this->field[iOwner];
                PhiHalo = this->field[iHalo];

                this->field[iNeighbour] = 0.5*(PhiOwner + PhiHalo);

            } else if (BCType == "empty") {

                this->field[iNeighbour] = this->field[iOwner];
            } else {

                printf("ERROR. No correct boundary condition type selected !!\n");
            }
        }
    }
}*/
