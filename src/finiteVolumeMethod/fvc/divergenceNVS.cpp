//
// Created by ruben on 31/01/24.
//

#include "fvc.h"


ScalarField fvc::divergenceNVS(const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs) {

    // Auxiliary variables
    int iOwner, iNeighbour, iPeriodicFace, iHalo;
    double gf;
    GeometricVector PhiF, PhiOwner, PhiNeighbour, PhiHalo, Sf, BCValue;
    std::string BCType;


    // Preallocate the divergence field
    ScalarField divergence;
    divergence.assign(theMesh.nInteriorElements, 0);


    // Loop over all the interior faces
    for (int i = 0; i < theMesh.nInteriorFaces; ++i) {

        gf = theMesh.faces[i].gf;
        Sf = theMesh.faces[i].Sf;

        iOwner = theMesh.faces[i].iOwner;
        iNeighbour = theMesh.faces[i].iNeighbour;

        PhiOwner = Phi[iOwner];
        PhiNeighbour = Phi[iNeighbour];

        PhiF = gf*PhiOwner + (1 - gf)*PhiNeighbour;

        divergence[iOwner] += PhiF*Sf;
        divergence[iNeighbour] -= PhiF*Sf;
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

                PhiF = {0,0,0};
            }

            divergence[iOwner] += PhiF*Sf;
        }
    }


    return divergence;
}