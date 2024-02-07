//
// Created by ruben on 6/02/24.
//

#include <cmath>
#include "initializeTurbulentChannelFlowReTau180.h"


VectorField initializeTurbulentChannelFlowReTau180(const PolyMesh& theMesh, double Lx, double Ly, double Lz, double nu) {

    // Auxiliary variables declaration
    double b = 5.5, kappa = 2.5, A = 0.1*(kappa*log(1/nu) + b);
    double uAvg, sinX, sinY, sinZ, cosX, cosY, cosZ;


    // Initialize the velocity field
    VectorField u;
    u.assign(theMesh.nInteriorElements, {0,0,0});


    // Loop for all the interior elements
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        double y = theMesh.elements[i].centroid.y - 1;

        // Average velocity computation
        if ((1 - fabs(y))/nu < 10) {

            uAvg = (1 - fabs(y))/nu;
        } else {

            uAvg = kappa*log((1 - fabs(y))/nu) + b;
        }

        // Perturbations computation
        sinX = sin(4*M_PI*theMesh.elements[i].centroid.x/Lx);
        sinY = sin(M_PI*y);
        sinZ = sin(2*M_PI*theMesh.elements[i].centroid.z/Lz);
        cosX = cos(4*M_PI*theMesh.elements[i].centroid.x/Lx);
        cosY = 1 + cos(2*M_PI*y/Ly);
        cosZ = cos(2*M_PI*theMesh.elements[i].centroid.z/Lz);

        // Velocity computation
        u[i].x = uAvg + A*cosX*sinY*sinZ*Lx/2;
        u[i].y = -A*sinX*cosY*sinZ*Ly/2;
        u[i].z = -A*sinX*sinY*cosZ*Lz/2;
    }


    return u;
}