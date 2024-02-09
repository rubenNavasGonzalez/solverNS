//
// Created by ruben on 8/02/24.
//

#include <cmath>
#include "modelWALE.h"
#include "../../../finiteVolumeMethod/fvc/fvc.h"


ScalarField modelWALE(const PolyMesh& theMesh, const VectorField& u, const VectorBoundaryConditions& uBCs) {

    // Auxiliary variables declaration
    TensorField gradU;
    double Delta, V2, QG, QS;
    double Cw = 0.569;


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
        V2 = gradU[i].invariantV2();
        QG = gradU[i].invariantQ();
        QS = gradU[i].symmetric().invariantQ();

        nut[i] = 2 * pow(Cw*Delta,2) * pow(V2/2 + 2./3*pow(QG,2),1.5) /
                (pow(-2*QS,2.5) + pow(V2/2 + 2./3*pow(QG,2),1.25));
    }


    return nut;
}