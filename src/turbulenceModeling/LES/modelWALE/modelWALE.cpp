//
// Created by ruben on 8/02/24.
//

#include <cmath>
#include "modelWALE.h"
#include "../../../finiteVolumeMethod/fvc/fvc.h"


ScalarField modelWALE(const PolyMesh& theMesh, const VectorField& u, const VectorBoundaryConditions& uBCs) {

    // Auxiliary variables declaration
    TensorField gradU;
    Tensor S, Omega, Sd;
    double Delta, QS, QSd;
    double Cw = sqrt(0.5);


    // Initialize the turbulent viscosity field
    ScalarField nut;
    nut.assign(theMesh.nInteriorElements, 0);


    // Compute the velocity gradient TensorField
    gradU = fvc::gradient(u, theMesh, uBCs);


    // Loop over all the interior elements
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {


        // Filter length computation
        Delta = pow(theMesh.elements[i].Vf, 1./3);


        // Compute the turbulent viscosity
        S = gradU[i].symmetric();
        Omega = gradU[i].skewSymmetric();
        Sd = S*S + Omega*Omega;
        QS = S.invariantQ();
        QSd = -0.5*(Sd*Sd).trace() + 1./6*pow(Sd.trace(), 2);

        nut[i] = pow(Cw*Delta, 2) * pow(-2*QSd, 1.5) / (pow(-2*QS, 2.5) + pow(-2*QSd, 1.25));
    }


    return nut;
}