//
// Created by ruben on 26/11/23.
//

#include "fvc.h"


VectorField fvc::convective(const ScalarField& mDot, const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs) {

    // Auxiliary variables
    int iOwner, iNeighbour, iPeriodicFace, iHalo;
    GeometricVector PhiF, PhiOwner, PhiNeighbour, PhiHalo, Sf, BCValue;
    std::string BCType;


    // Preallocate the divergence field
    VectorField convective;
    convective.assign(theMesh.nInteriorElements, {0,0,0});


    // Loop over all the interior faces
    for (int i = 0; i < theMesh.nInteriorFaces; ++i) {

        // Get the owner and neighbour element indices
        iOwner = theMesh.faces[i].iOwner;
        iNeighbour = theMesh.faces[i].iNeighbour;

        // Get the owner and neighbour Phi value
        PhiOwner = Phi[iOwner];
        PhiNeighbour = Phi[iNeighbour];

        // Interpolate the Phi value at the element face (a symmetry-preserving scheme is used)
        PhiF = 0.5*(PhiOwner + PhiNeighbour);

        // Compute the convective term
        convective[iOwner] += mDot[i]*PhiF;
        convective[iNeighbour] -= mDot[i]*PhiF;
    }


    // Loop over all the boundaries
    for (int k = 0; k < theMesh.nBoundaries; ++k) {

        BCType = PhiBCs.type[k];
        BCValue = PhiBCs.value[k];

        // Loop over all the boundary faces of the k boundary
        for (int i = theMesh.boundaries[k].startFace; i < theMesh.boundaries[k].startFace + theMesh.boundaries[k].nBoundaryFaces; ++i) {

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

            convective[iOwner] += mDot[i]*PhiF;
        }
    }


    // Loop over all the interior elements
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        convective[i] /= theMesh.elements[i].Vf;
    }


    return convective;
}


ScalarField fvc::convective(const ScalarField& mDot, const ScalarField& Phi, const PolyMesh& theMesh, const ScalarBoundaryConditions& PhiBCs) {

    // Auxiliary variables
    int iOwner, iNeighbour, iPeriodicFace, iHalo;
    double PhiF, PhiOwner, PhiNeighbour, PhiHalo, BCValue;
    GeometricVector Sf;
    std::string BCType;


    // Preallocate the divergence field
    ScalarField convective;
    convective.assign(theMesh.nInteriorElements, 0);


    // Loop over all the interior faces
    for (int i = 0; i < theMesh.nInteriorFaces; ++i) {

        // Get the owner and neighbour element indices
        iOwner = theMesh.faces[i].iOwner;
        iNeighbour = theMesh.faces[i].iNeighbour;

        // Get the owner and neighbour Phi value
        PhiOwner = Phi[iOwner];
        PhiNeighbour = Phi[iNeighbour];

        // Interpolate the Phi value at the element face (a symmetry-preserving scheme is used)
        PhiF = 0.5*(PhiOwner + PhiNeighbour);

        // Compute the convective term
        convective[iOwner] += mDot[i]*PhiF;
        convective[iNeighbour] -= mDot[i]*PhiF;
    }


    // Loop over all the boundaries
    for (int k = 0; k < theMesh.nBoundaries; ++k) {

        BCType = PhiBCs.type[k];
        BCValue = PhiBCs.value[k];

        // Loop over all the boundary faces of the k boundary
        for (int i = theMesh.boundaries[k].startFace; i < theMesh.boundaries[k].startFace + theMesh.boundaries[k].nBoundaryFaces; ++i) {

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

            convective[iOwner] += mDot[i]*PhiF;
        }
    }


    // Loop over all the interior elements
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        convective[i] /= theMesh.elements[i].Vf;
    }


    return convective;
}