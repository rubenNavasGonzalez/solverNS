//
// Created by ruben on 6/02/24.
//

#include "VectorField.h"
#include "initializeTurbulentChannelFlowReTau180.h"
#include "initializeFromSelectedTime.h"
#include "interpolateFromDifferentMesh.h"


void VectorField::initialize(const PolyMesh& theMesh, double t, double Lx, double Ly, double Lz,  double nu, int mode) {

    if (t == 0) {

        // Initialize the velocity field according to the selected case
        *this = initializeTurbulentChannelFlowReTau180(theMesh, Lx, Ly, Lz, nu);

    } else {

        if (mode == 0) {

            // Initialize the velocity field from time t and equal mesh
            *this = initializeFromSelectedTime(theMesh, t);

        } else if (mode == 1) {

            // Initialize the velocity field from time t and different mesh
            *this = interpolateFromDifferentMesh(theMesh, t);

        } else {

            printf("Error. No correct initialization mode selected. Valid ones are 0 (start from time t) and 1 "
                   "(interpolate from different mesh at time t). \n");
            std::exit(EXIT_FAILURE);
        }
    }
}