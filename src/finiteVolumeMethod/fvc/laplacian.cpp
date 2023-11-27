//
// Created by ruben on 26/11/23.
//

#include "laplacian.h"


VectorField fvc::laplacian(const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs) {

    // Auxiliary variables
    int iOwner, iNeighbour, iPeriodicFace;
    double gf, DivPhiF, DivPhiOwner, DivPhiNeighbour;
    GeometricVector Sf, BCValue;
    std::string BCType;


    // Preallocate the laplacian field
    VectorField laplacian;
    laplacian.initialize(theMesh.nElements);


    // Compute the divergence field
    ScalarField DivPhi = fvc::divergence(Phi, theMesh, PhiBCs);


    // Loop over all the interior faces
    for (int i = 0; i < theMesh.nInteriorFaces; ++i) {

        gf = theMesh.faces[i].gf;
        Sf = theMesh.faces[i].Sf;

        iOwner = theMesh.faces[i].iOwner;
        iNeighbour = theMesh.faces[i].iNeighbour;

        DivPhiOwner = DivPhi.field[iOwner];
        DivPhiNeighbour = DivPhi.field[iNeighbour];

        DivPhiF = gf*DivPhiOwner + (1 - gf)*DivPhiNeighbour;

        laplacian.field[iOwner] += DivPhiF*Sf;
        laplacian.field[iNeighbour] -= DivPhiF*Sf;
    }


    // Loop over all the boundaries
    for (int k = 0; k < theMesh.nBoundaries; ++k) {

        BCType = PhiBCs.type[k];

        // Loop over all the boundary faces of the k boundary
        for (int i = theMesh.boundaries[k].startFace; i < theMesh.boundaries[k].startFace + theMesh.boundaries[k].nBoundaryFaces; ++i) {

            Sf = theMesh.faces[i].Sf;

            iOwner = theMesh.faces[i].iOwner;
            iNeighbour = theMesh.faces[i].iNeighbour;
            iPeriodicFace = theMesh.faces[i].iPeriodicFace;

            if (BCType == "zeroGradient") {

                DivPhiF = 0;
            } else if (BCType == "periodic") {

                DivPhiOwner = DivPhi.field[iNeighbour];
                DivPhiNeighbour = DivPhi.field[theMesh.faces[iPeriodicFace].iNeighbour];

                DivPhiF = 0.5*(DivPhiOwner + DivPhiNeighbour);
            } else if (BCType == "empty") {

                DivPhiF = 0;
            }

            laplacian.field[iOwner] += DivPhiF*Sf;
        }
    }


    // Loop over all the interior elements
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        laplacian.field[i] /= theMesh.elements[i].Vf;
    }


    // Loop over all the boundary faces
    for (int i = theMesh.nInteriorFaces; i < theMesh.nFaces; ++i) {

        iOwner = theMesh.faces[i].iOwner;
        iNeighbour = theMesh.faces[i].iNeighbour;

        laplacian.field[iNeighbour] = laplacian.field[iOwner];
    }


    return laplacian;
}