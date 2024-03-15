//
// Created by ruben on 3/03/24.
//

#include <cmath>
#include "initializeTaylorGreenVortex.h"


VectorField initializeTaylorGreenVortex(const PolyMesh& theMesh, double Lx, double Ly, double Lz) {


    // Auxiliary variables definition
    double x, y, z;


    // Initialize the velocity field
    VectorField u;
    u.assign(theMesh.nInteriorElements, {0,0,0});


    // Loop over all the interior elements
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {


        // Get the coordinates of the centroid
        x = theMesh.elements[i].centroid.x;
        y = theMesh.elements[i].centroid.y;
        z = theMesh.elements[i].centroid.z;


        // Assign the corresponding velocity field values
        u[i].x = sin(x)*cos(y)*cos(z);
        u[i].y = -cos(x)*sin(y)*cos(z);
        u[i].z = 0;
    }


    return u;
}
