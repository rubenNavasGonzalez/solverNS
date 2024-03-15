//
// Created by ruben on 8/02/24.
//

#include <cmath>
#include "modelSmagorinsky.h"
#include "../../../finiteVolumeMethod/fvc/fvc.h"


ScalarField modelSmagorinsky(const PolyMesh& theMesh, const VectorField& u, const VectorBoundaryConditions& uBCs) {

    // Auxiliary variables declaration
    TensorField gradU;
    double Delta, QS;
    double Cs = 0.165;


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
        QS = gradU[i].symmetric().invariantQ();
        nut[i] = 2 * pow(Cs*Delta, 2) * sqrt(-QS);
    }


    return nut;
}