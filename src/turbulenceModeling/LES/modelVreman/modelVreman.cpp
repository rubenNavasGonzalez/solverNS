//
// Created by ruben on 8/02/24.
//

#include <cmath>
#include "modelVreman.h"
#include "../../../finiteVolumeMethod/fvc/fvc.h"


ScalarField modelVreman(const PolyMesh& theMesh, const VectorField& u, const VectorBoundaryConditions& uBCs) {

    // Auxiliary variables declaration
    TensorField gradU;
    Tensor GGT;
    double Delta;
    double Cvr = sqrt(0.07);


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
        GGT = gradU[i] * (gradU[i].transpose());
        nut[i] = 2 * pow(Cvr*Delta,2) * pow( GGT.invariantQ() / GGT.invariantP() , 0.5);
    }


    return nut;
}