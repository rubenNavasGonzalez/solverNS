//
// Created by ruben on 11/12/23.
//

#include "RhieChowInterpolation.h"
#include "../../finiteVolumeMethod/fvc/fvc.h"


ScalarField RhieChowInterpolation(const VectorField& u, const ScalarField& p, const PolyMesh& theMesh, double DeltaT, const VectorBoundaryConditions& uBCs, const ScalarBoundaryConditions& pBCs) {

    // Auxiliary variables definition
    int iOwner, iNeighbour, iPeriodicFace, iHalo;
    double gf, dONMag, mDotAvg, mDotCorr;
    GeometricVector Sf, uF, gradPF, BCValue;
    std::string BCType;


    // Preallocate mass flow field
    ScalarField mDot;
    mDot.assign(theMesh.nFaces, 0);


    // Compute the pressure gradient
    VectorField gradP = fvc::gradient(p, theMesh, pBCs);


    // Loop over all the interior faces
    for (int i = 0; i < theMesh.nInteriorFaces; ++i) {

        iOwner = theMesh.faces[i].iOwner;
        iNeighbour = theMesh.faces[i].iNeighbour;

        gf = theMesh.faces[i].gf;
        Sf = theMesh.faces[i].Sf;
        dONMag = theMesh.faces[i].dONMag;

        uF = gf*u[iOwner] + (1 - gf)*u[iNeighbour];
        mDotAvg = 1*uF*Sf;

        gradPF = gf*gradP[iOwner] + (1 - gf)*gradP[iNeighbour];
        mDotCorr = DeltaT*(gradPF*Sf - (p[iNeighbour] - p[iOwner])/dONMag*Sf.mag());

        mDot[i] = mDotAvg + mDotCorr;
    }


    // Loop over all the boundaries
    for (int k = 0; k < theMesh.nBoundaries; ++k) {

        BCType = uBCs.type[k];
        BCValue = uBCs.value[k];

        // Loop over all the boundary faces of the k boundary
        for (int i = theMesh.boundaries[k].startFace; i < theMesh.boundaries[k].startFace + theMesh.boundaries[k].nBoundaryFaces; ++i) {

            Sf = theMesh.faces[i].Sf;
            dONMag = theMesh.faces[i].dONMag;
            iOwner = theMesh.faces[i].iOwner;

            if (BCType == "fixedValue") {

                mDot[i] = 1*BCValue*Sf;
            } else if (BCType == "zeroGradient") {

                mDot[i] = 1*u[iOwner]*Sf;
            } else if (BCType == "periodic") {

                iPeriodicFace = theMesh.faces[i].iPeriodicFace;
                iHalo = theMesh.faces[iPeriodicFace].iOwner;

                uF = 0.5*(u[iOwner] + u[iHalo]);
                mDotAvg = 1*uF*Sf;

                gradPF = 0.5*(gradP[iOwner] + gradP[iHalo]);
                mDotCorr = DeltaT*(gradPF*Sf - (p[iHalo] - p[iOwner])/(2*dONMag)*Sf.mag());

                mDot[i] = mDotAvg + mDotCorr;

            } else {

                mDot[i] = 0;
            }
        }
    }


    return mDot;
}