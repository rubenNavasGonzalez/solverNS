//
// Created by ruben on 27/11/23.
//

#include "laplacianOrthogonal.h"


VectorField fvc::laplacianOrthogonal(const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs) {

    // Auxiliary variables
    double dONMag, SfMag;
    int iOwner, iNeighbour, iHalo, iPeriodicFace;
    GeometricVector PhiOwner, PhiNeighbour, PhiHalo, DivPhiF, BCValue;
    std::string BCType;


    // Preallocate the laplacian field
    VectorField laplacian;
    laplacian.initialize(theMesh.nElements);


    // Loop over all the interior faces
    for (int i = 0; i < theMesh.nInteriorFaces; ++i) {

        dONMag = theMesh.faces[i].dONMag;
        SfMag = theMesh.faces[i].SfMag;

        iOwner = theMesh.faces[i].iOwner;
        iNeighbour = theMesh.faces[i].iNeighbour;

        PhiOwner = Phi.field[iOwner];
        PhiNeighbour = Phi.field[iNeighbour];

        DivPhiF = (PhiNeighbour - PhiOwner)/dONMag;

        laplacian.field[iOwner] += DivPhiF*SfMag;
        laplacian.field[iNeighbour] -= DivPhiF*SfMag;
    }


    // Loop over all the boundaries
    for (int k = 0; k < theMesh.nBoundaries; ++k) {

        BCType = PhiBCs.type[k];
        BCValue = PhiBCs.value[k];


        // Loop over all the boundary faces of the k_th boundary
        for (int i = theMesh.boundaries[k].startFace; i < theMesh.boundaries[k].startFace + theMesh.boundaries[k].nBoundaryFaces; ++i) {


            if (BCType == "fixedValue") {

                DivPhiF = BCValue;
            } else if (BCType == "zeroGradient" || BCType == "empty") {

                DivPhiF = {0,0,0};
            } else if (BCType == "periodic") {

                dONMag = theMesh.faces[i].dONMag;
                SfMag = theMesh.faces[i].SfMag;

                iPeriodicFace = theMesh.faces[i].iPeriodicFace;
                iOwner = theMesh.faces[i].iOwner;
                iHalo = theMesh.faces[iPeriodicFace].iOwner;

                PhiOwner = Phi.field[iOwner];
                PhiHalo = Phi.field[iHalo];

                DivPhiF = (PhiHalo - PhiOwner)/(2*dONMag);

                } else {

                printf("ERROR. No correct boundary condition type selected !!\n");
            }

            laplacian.field[iOwner] += DivPhiF*SfMag;
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