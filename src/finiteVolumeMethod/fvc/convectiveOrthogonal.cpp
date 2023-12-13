//
// Created by ruben on 26/11/23.
//

#include "fvc.h"
#include "../../interpolation/convectiveSchemes/convectiveSchemesOrthogonal.h"


VectorField fvc::convectiveOrthogonal(const ScalarField& mDot, const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs, const std:: string& scheme) {

    // Auxiliary variables
    int iOwner, iNeighbour, iOwnerFar, iNeighbourFar, iPeriodicFace, iHalo;
    Node pOwner, pNeighbour, pOwnerFar, pNeighbourFar, pF;
    GeometricVector PhiF, PhiOwner, PhiNeighbour, PhiOwnerFar, PhiNeighbourFar, PhiHalo, Sf, BCValue;
    std::string BCType;


    // Preallocate the divergence field
    VectorField convective;
    convective.initialize(theMesh.nElements);


    // Loop over all the interior faces
    for (int i = 0; i < theMesh.nInteriorFaces; ++i) {

        Sf = theMesh.faces[i].Sf;

        iOwner = theMesh.faces[i].iOwner;
        iNeighbour = theMesh.faces[i].iNeighbour;
        iOwnerFar = theMesh.faces[i].iOwnerFar;
        iNeighbourFar = theMesh.faces[i].iNeighbourFar;

        pOwner = theMesh.elements[iOwner].centroid;
        pNeighbour = theMesh.elements[iNeighbour].centroid;
        pOwnerFar = theMesh.elements[iOwnerFar].centroid;
        pNeighbourFar = theMesh.elements[iNeighbourFar].centroid;
        pF = theMesh.faces[i].centroid;

        PhiOwner = Phi[iOwner];
        PhiNeighbour = Phi[iNeighbour];
        PhiOwnerFar = Phi[iOwnerFar];
        PhiNeighbourFar = Phi[iNeighbourFar];

        PhiF = convectiveSchemesOrthogonal(PhiOwner, pOwner, PhiOwnerFar, pOwnerFar, PhiNeighbour, pNeighbour, PhiNeighbourFar,
                                            pNeighbourFar, pF, mDot[i], Sf, scheme);

        convective.field[iOwner] += mDot[i]*PhiF;
        convective.field[iNeighbour] -= mDot[i]*PhiF;
    }


    // Loop over all the boundaries
    for (int k = 0; k < theMesh.nBoundaries; ++k) {

        BCType = PhiBCs.type[k];
        BCValue = PhiBCs.value[k];

        // Loop over all the boundary faces of the k boundary
        for (int i = theMesh.boundaries[k].startFace; i < theMesh.boundaries[k].startFace + theMesh.boundaries[k].nBoundaryFaces; ++i) {

            if (BCType == "fixedValue") {

                PhiF = BCValue;
            } else if (BCType == "zeroGradient") {

                PhiF = Phi[iOwner];
            } else if (BCType == "periodic") {

                iOwner = theMesh.faces[i].iOwner;
                iPeriodicFace = theMesh.faces[i].iPeriodicFace;
                iHalo = theMesh.faces[iPeriodicFace].iOwner;

                PhiOwner = Phi[iOwner];
                PhiHalo = Phi[iHalo];

                PhiF = 0.5*(PhiOwner + PhiHalo);
            } else {

                PhiF = {0,0,0};
            }

            convective.field[iOwner] += mDot[i]*PhiF;
        }
    }


    // Loop over all the interior elements
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        convective.field[i] /= theMesh.elements[i].Vf;
    }


    // Loop over all the boundary faces
    for (int i = theMesh.nInteriorFaces; i < theMesh.nFaces; ++i) {

        iOwner = theMesh.faces[i].iOwner;
        iNeighbour = theMesh.faces[i].iNeighbour;

        convective.field[iNeighbour] = convective[iOwner];
    }


    return convective;
}