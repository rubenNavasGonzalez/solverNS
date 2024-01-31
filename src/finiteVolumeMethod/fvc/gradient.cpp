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
    gradient.assign(theMesh.nInteriorElements, {0,0,0});


    // Loop over all the interior faces
    for (int i = 0; i < theMesh.nInteriorFaces; ++i) {

        gf = theMesh.faces[i].gf;
        Sf = theMesh.faces[i].Sf;

        iOwner = theMesh.faces[i].iOwner;
        iNeighbour = theMesh.faces[i].iNeighbour;

        PhiOwner = Phi[iOwner];
        PhiNeighbour = Phi[iNeighbour];

        PhiF = gf*PhiOwner + (1 - gf)*PhiNeighbour;

        gradient[iOwner] += PhiF*Sf;
        gradient[iNeighbour] -= PhiF*Sf;
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

            gradient[iOwner] += PhiF*Sf;
        }
    }


    // Loop over all the interior elements
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        gradient[i] /= theMesh.elements[i].Vf;
    }


    return gradient;
}



TensorField fvc::gradient(const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs) {

    // Auxiliary variables declarations
    ScalarField PhiX, PhiY, PhiZ;
    VectorField gradX, gradY, gradZ;
    ScalarBoundaryConditions PhiXBCs, PhiYBCs, PhiZBCs;


    // Preallocate the gradient field
    TensorField gradient;
    gradient.resize(theMesh.nInteriorElements);


    // Split the velocity vector field in the gradient field of its components
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        PhiX.push_back(Phi[i].x);
        PhiY.push_back(Phi[i].y);
        PhiZ.push_back(Phi[i].z);
    }


    // Split the velocity boundary conditions in scalar boundary conditions for each component
    for (int i = 0; i < PhiBCs.type.size(); ++i) {

        PhiXBCs.type.push_back(PhiBCs.type[i]);
        PhiXBCs.value.push_back(PhiBCs.value[i].x);

        PhiYBCs.type.push_back(PhiBCs.type[i]);
        PhiYBCs.value.push_back(PhiBCs.value[i].y);

        PhiZBCs.type.push_back(PhiBCs.type[i]);
        PhiZBCs.value.push_back(PhiBCs.value[i].z);
    }


    // Compute the gradient of each velocity component field
    gradX = fvc::gradient(PhiX, theMesh, PhiXBCs);
    gradY = fvc::gradient(PhiY, theMesh, PhiYBCs);
    gradZ = fvc::gradient(PhiZ, theMesh, PhiZBCs);


    // Assemble the gradient tensor field
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        gradient[i][0][0] = gradX[i].x;
        gradient[i][0][1] = gradY[i].x;
        gradient[i][0][2] = gradZ[i].x;

        gradient[i][1][0] = gradX[i].y;
        gradient[i][1][1] = gradY[i].y;
        gradient[i][1][2] = gradZ[i].y;

        gradient[i][2][0] = gradX[i].z;
        gradient[i][2][1] = gradY[i].z;
        gradient[i][2][2] = gradZ[i].z;
    }


    return gradient;
}