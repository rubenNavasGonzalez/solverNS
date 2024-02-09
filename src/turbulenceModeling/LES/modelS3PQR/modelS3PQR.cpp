//
// Created by ruben on 9/02/24.
//

#include <cmath>
#include "modelS3PQR.h"
#include "../../../finiteVolumeMethod/fvc/fvc.h"


ScalarField modelS3PQR(const PolyMesh& theMesh, const VectorField& u, const VectorBoundaryConditions& uBCs, modelLES turbulenceModel) {

    // Auxiliary variables declaration
    TensorField gradU;
    double Delta, C_s3pqr, p;
    Tensor GGT;
    double P_GGT, Q_GGT, R_GGT;


    // Assign constant model and exponent degree as a function of the turbulence model
    switch (turbulenceModel) {

        case S3PQ: {

            C_s3pqr = 0.458;
            p = -2.5;
        }
        case S3PR: {

            C_s3pqr = 0.458;
            p = -1;
        }
        case S3QR: {

            C_s3pqr = 0.458;
            p = 0;
        }
    }


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
        GGT = gradU[i] * gradU[i].transpose();
        P_GGT = GGT.invariantP();
        Q_GGT = GGT.invariantQ();
        R_GGT = GGT.invariantR();

        nut[i] = 2 * pow(C_s3pqr*Delta,2) * pow(P_GGT,p) * pow(Q_GGT,-(p + 1)) * pow(R_GGT,(p + 2.5)/3);
    }


    return nut;
}