//
// Created by ruben on 7/12/23.
//

#include "fvc.h"


VectorField fvc::gradient(const ScalarField& Phi, const PolyMesh& theMesh, const ScalarBoundaryConditions& PhiBCs) {

    // Auxiliary variables
    int iOwner, iNeighbour, iPeriodicFace, iHalo;
    double gf, PhiF, PhiOwner, PhiNeighbour, PhiHalo, BCValue;
    GeometricVector Sf;
    std::string BCType;


    // Preallocate the gradient field
    VectorField gradient;
    //gradient.initialize(theMesh.nElements);
    gradient.initialize(theMesh.nInteriorElements);


    // Loop over all the interior faces
    for (int i = 0; i < theMesh.nInteriorFaces; ++i) {

        gf = theMesh.faces[i].gf;
        Sf = theMesh.faces[i].Sf;

        iOwner = theMesh.faces[i].iOwner;
        iNeighbour = theMesh.faces[i].iNeighbour;

        PhiOwner = Phi[iOwner];
        PhiNeighbour = Phi[iNeighbour];

        PhiF = gf*PhiOwner + (1 - gf)*PhiNeighbour;

        gradient.field[iOwner] += PhiF*Sf;
        gradient.field[iNeighbour] -= PhiF*Sf;
    }


    // Loop over all the boundaries
    for (int k = 0; k < theMesh.nBoundaries; ++k) {

        BCType = PhiBCs.type[k];
        BCValue = PhiBCs.value[k];

        // Loop over all the boundary faces of the k boundary
        for (int i = theMesh.boundaries[k].startFace; i < theMesh.boundaries[k].startFace + theMesh.boundaries[k].nBoundaryFaces; ++i) {

            Sf = theMesh.faces[i].Sf;
            iOwner = theMesh.faces[i].iOwner;

            if (BCType == "fixedValue") {

                PhiF = BCValue;
            } else if (BCType == "zeroGradient") {

                PhiF = Phi[iOwner];
            } else if (BCType == "periodic") {

                iPeriodicFace = theMesh.faces[i].iPeriodicFace;
                iHalo = theMesh.faces[iPeriodicFace].iOwner;

                PhiOwner = Phi[iOwner];
                PhiHalo = Phi[iHalo];

                PhiF = 0.5*(PhiOwner + PhiHalo);

            } else {

                PhiF = 0;
            }

            gradient.field[iOwner] += PhiF*Sf;
        }
    }


    // Loop over all the interior elements
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        gradient.field[i] /= theMesh.elements[i].Vf;
    }


    // Loop over all the boundary faces
    /*for (int i = theMesh.nInteriorFaces; i < theMesh.nFaces; ++i) {

        iOwner = theMesh.faces[i].iOwner;
        iNeighbour = theMesh.faces[i].iNeighbour;

        gradient.field[iNeighbour] = gradient[iOwner];
    }*/


    return gradient;
}