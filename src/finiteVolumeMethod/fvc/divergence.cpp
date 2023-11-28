//
// Created by ruben on 26/11/23.
//

#include "divergence.h"


ScalarField fvc::divergence(const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs) {

    // Auxiliary variables
    int iOwner, iNeighbour, iPeriodicFace;
    double gf;
    GeometricVector PhiF, PhiOwner, PhiNeighbour, Sf, BCValue;
    std::string BCType;


    // Preallocate the divergence field
    ScalarField divergence;
    divergence.initialize(theMesh.nElements);


    // Loop over all the interior faces
    for (int i = 0; i < theMesh.nInteriorFaces; ++i) {

        gf = theMesh.faces[i].gf;
        Sf = theMesh.faces[i].Sf;

        iOwner = theMesh.faces[i].iOwner;
        iNeighbour = theMesh.faces[i].iNeighbour;

        PhiOwner = Phi.field[iOwner];
        PhiNeighbour = Phi.field[iNeighbour];

        PhiF = gf*PhiOwner + (1 - gf)*PhiNeighbour;

        divergence.field[iOwner] += PhiF*Sf;
        divergence.field[iNeighbour] -= PhiF*Sf;
    }


    // Loop over all the boundaries
    for (int k = 0; k < theMesh.nBoundaries; ++k) {

        BCType = PhiBCs.type[k];
        BCValue = PhiBCs.value[k];

        // Loop over all the boundary faces of the k boundary
        for (int i = theMesh.boundaries[k].startFace; i < theMesh.boundaries[k].startFace + theMesh.boundaries[k].nBoundaryFaces; ++i) {

            Sf = theMesh.faces[i].Sf;

            iOwner = theMesh.faces[i].iOwner;
            iNeighbour = theMesh.faces[i].iNeighbour;
            iPeriodicFace = theMesh.faces[i].iPeriodicFace;

            if (BCType == "fixedValue") {

                PhiF = BCValue;
            } else if (BCType == "zeroGradient") {

                PhiF = Phi.field[iNeighbour];
            } else if (BCType == "periodic") {

                PhiOwner = Phi.field[iNeighbour];
                PhiNeighbour = Phi.field[theMesh.faces[iPeriodicFace].iNeighbour];

                PhiF = 0.5*(PhiOwner + PhiNeighbour);
            } else {

                PhiF = {0,0,0};
            }

            divergence.field[iOwner] += PhiF*Sf;
        }
    }


    // Loop over all the interior elements
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        divergence.field[i] /= theMesh.elements[i].Vf;
    }


    // Loop over all the boundary faces
    for (int i = theMesh.nInteriorFaces; i < theMesh.nFaces; ++i) {

        iOwner = theMesh.faces[i].iOwner;
        iNeighbour = theMesh.faces[i].iNeighbour;

        divergence.field[iNeighbour] = divergence.field[iOwner];
    }

    return divergence;
}