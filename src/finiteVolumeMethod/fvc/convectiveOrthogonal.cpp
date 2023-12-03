//
// Created by ruben on 26/11/23.
//

#include "convectiveOrthogonal.h"
#include "../../interpolation/convectiveSchemes/convectiveSchemesOrthogonal.h"


VectorField fvc::convectiveOrthogonal(const ScalarField& mDot, const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs, const std:: string& scheme) {

    // Auxiliary variables
    int iOwner, iNeighbour, iOwnerFar, iNeighbourFar, iPeriodicFace, iHalo;
    double gf, mDotVal, pOwner, pNeighbour, pOwnerFar, pNeighbourFar, pFace;
    GeometricVector PhiF, PhiFAvg, PhiOwner, PhiNeighbour, PhiOwnerFar, PhiNeighbourFar, PhiHalo, Sf, BCValue;
    std::string BCType;


    // Preallocate the divergence field
    VectorField convective;
    convective.initialize(theMesh.nElements);


    // Loop over all the interior faces
    for (int i = 0; i < theMesh.nInteriorFaces; ++i) {

        gf = theMesh.faces[i].gf;
        Sf = theMesh.faces[i].Sf;

        iOwner = theMesh.faces[i].iOwner;
        iNeighbour = theMesh.faces[i].iNeighbour;
        iOwnerFar = theMesh.faces[i].iOwnerFar;
        iNeighbourFar = theMesh.faces[i].iNeighbourFar;

        PhiOwner = Phi.field[iOwner];
        PhiNeighbour = Phi.field[iNeighbour];
        PhiOwnerFar = Phi.field[iOwnerFar];
        PhiNeighbourFar = Phi.field[iNeighbourFar];

        PhiF = gf*PhiOwner + (1 - gf)*PhiNeighbour;
        mDotVal = 1*PhiF*Sf;

        if (Sf.x != 0) {

            pOwner = theMesh.elements[iOwner].centroid.x;
            pNeighbour = theMesh.elements[iNeighbour].centroid.x;
            pOwnerFar = theMesh.elements[iOwnerFar].centroid.x;
            pNeighbourFar = theMesh.elements[iNeighbourFar].centroid.x;
            pFace = theMesh.faces[i].centroid.x;

            PhiF = convectiveSchemesOrthogonal(PhiOwner, pOwner, PhiOwnerFar, pOwnerFar, PhiNeighbour,
                                               pNeighbour, PhiNeighbourFar, pNeighbourFar, pFace, mDotVal, scheme);

        } else if (Sf.y != 0) {

            pOwner = theMesh.elements[iOwner].centroid.y;
            pNeighbour = theMesh.elements[iNeighbour].centroid.y;
            pOwnerFar = theMesh.elements[iOwnerFar].centroid.y;
            pNeighbourFar = theMesh.elements[iNeighbourFar].centroid.y;
            pFace = theMesh.faces[i].centroid.y;

            PhiF = convectiveSchemesOrthogonal(PhiOwner, pOwner, PhiOwnerFar, pOwnerFar, PhiNeighbour,
                                               pNeighbour, PhiNeighbourFar, pNeighbourFar, pFace, mDotVal, scheme);

        } else {

            pOwner = theMesh.elements[iOwner].centroid.z;
            pNeighbour = theMesh.elements[iNeighbour].centroid.z;
            pOwnerFar = theMesh.elements[iOwnerFar].centroid.z;
            pNeighbourFar = theMesh.elements[iNeighbourFar].centroid.z;
            pFace = theMesh.faces[i].centroid.z;

            PhiF = convectiveSchemesOrthogonal(PhiOwner, pOwner, PhiOwnerFar, pOwnerFar, PhiNeighbour,
                                               pNeighbour, PhiNeighbourFar, pNeighbourFar, pFace, mDotVal, scheme);

        }

        convective.field[iOwner] += mDotVal*PhiF;
        convective.field[iNeighbour] -= mDotVal*PhiF;
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
            iHalo = theMesh.faces[iPeriodicFace].iOwner;

            if (BCType == "fixedValue") {

                PhiF = BCValue;
                mDotVal = 1*PhiF*Sf;
            } else if (BCType == "zeroGradient") {

                PhiF = Phi.field[iOwner];
                mDotVal = 1*PhiF*Sf;
            } else if (BCType == "periodic") {

                PhiOwner = Phi.field[iOwner];
                PhiHalo = Phi.field[iHalo];

                PhiF = 0.5*(PhiOwner + PhiHalo);
                mDotVal = 1*PhiF*Sf;

                if (mDotVal >= 0) {

                    PhiF = PhiOwner;
                } else {

                    PhiF = PhiHalo;
                }
            } else {

                PhiF = {0,0,0};
                mDotVal = 0;
            }

            /*if (BCType != "empty") {

                PhiF = Phi.field[iNeighbour];
                mDotVal = 1*PhiF*Sf;
            } else {

                PhiF = {0,0,0};
                mDotVal = 0;
            }*/

            convective.field[iOwner] += mDotVal*PhiF;
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

        convective.field[iNeighbour] = convective.field[iOwner];
    }


    return convective;
}