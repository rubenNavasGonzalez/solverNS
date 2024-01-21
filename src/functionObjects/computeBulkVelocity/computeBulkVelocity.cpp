//
// Created by ruben on 26/12/23.
//

#include "computeBulkVelocity.h"


double computeBulkVelocity(const VectorField& u, const PolyMesh& theMesh, const VectorBoundaryConditions& uBCs, int k) {

    // Auxiliary variables initialization
    int iOwner, iPeriodicFace, iHalo;
    GeometricVector Sf, uF;
    std::string BCType = uBCs.type[k];
    GeometricVector BCValue = uBCs.value[k];


    // Mass flow and surface initialization
    double mDot{}, SfBoundary{};


    // Loop over all the boundary faces of the k boundary
    for (int i = theMesh.boundaries[k].startFace; i < theMesh.boundaries[k].startFace + theMesh.boundaries[k].nBoundaryFaces; ++i) {

        Sf = theMesh.faces[i].Sf;

        iOwner = theMesh.faces[i].iOwner;
        iPeriodicFace = theMesh.faces[i].iPeriodicFace;
        iHalo = theMesh.faces[iPeriodicFace].iOwner;

        if (BCType == "fixedValue") {

            mDot += 1*BCValue*Sf;
            SfBoundary += Sf.mag();

        } else if (BCType == "zeroGradient") {

            mDot += 1*u[iOwner]*Sf;
            SfBoundary += Sf.mag();

        } else if (BCType == "periodic") {

            uF = 0.5*(u[iOwner] + u[iHalo]);
            mDot += 1*uF*Sf;
            SfBoundary += Sf.mag();

        } else {

            SfBoundary += Sf.mag();
        }
    }


    return mDot/SfBoundary;
}