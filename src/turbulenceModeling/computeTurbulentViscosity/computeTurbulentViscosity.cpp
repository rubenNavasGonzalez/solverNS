//
// Created by ruben on 8/02/24.
//

#include <string>
#include "computeTurbulentViscosity.h"


ScalarField computeTurbulentViscosity(const PolyMesh& theMesh, const VectorField& u, const VectorBoundaryConditions& uBCs, modelLES turbulenceModel) {

    // Assign the turbulence viscosity depending on the turbulence model selected
    switch (turbulenceModel) {

        case None: {

            ScalarField nut;
            nut.assign(theMesh.nInteriorElements, 0);

            return nut;
        }
        case Smagorinsky: {
            return modelSmagorinsky(theMesh, u, uBCs);
        }
        case WALE: {
            return modelWALE(theMesh, u, uBCs);
        }
        case Vreman: {
            return modelVreman(theMesh, u, uBCs);
        }
        case Verstappen: {
            return modelVerstappen(theMesh, u, uBCs);
        }
        case S3PQ: {
            return modelS3PQR(theMesh, u, uBCs, turbulenceModel);
        }
        case S3PR: {
            return modelS3PQR(theMesh, u, uBCs, turbulenceModel);
        }
        case S3QR: {
            return modelS3PQR(theMesh, u, uBCs, turbulenceModel);
        }
        default: {

            printf("Error. Invalid turbulence model selected. \n");
            std::exit(EXIT_FAILURE);
        }
    }
}