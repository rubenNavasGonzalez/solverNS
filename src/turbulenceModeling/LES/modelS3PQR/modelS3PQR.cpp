//
// Created by ruben on 9/02/24.
//

#include <cmath>
#include "modelS3PQR.h"
#include "../../../finiteVolumeMethod/fvc/fvc.h"


ScalarField modelS3PQR(const PolyMesh& theMesh, const VectorField& u, const VectorBoundaryConditions& uBCs, modelLES turbulenceModel) {

    // Auxiliary variables declaration
    TensorField gradU;
    double Delta, p, P_GGT, Q_GGT, R_GGT;
    double C_s3pqr = sqrt(3) * sqrt(0.07);
    Tensor GGT;


    // Assign the and exponent degree depending on the turbulence model
    switch (turbulenceModel) {

        case S3PQ: {

            p = -2.5;
            //C_s3pqr = 0.572;
            break;
        }
        case S3PR: {

            p = -1;
            //C_s3pqr = 0.709;
            break;
        }
        case S3QR: {

            p = 0;
            //C_s3pqr = 0.762;
            break;
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