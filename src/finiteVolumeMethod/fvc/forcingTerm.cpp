//
// Created by ruben on 26/12/23.
//

#include <cmath>
#include "fvc.h"


VectorField fvc::forcingTerm(const GeometricVector& Fe, const PolyMesh& theMesh) {

    // Forcing term vector initialization
    VectorField F;
    F.assign(theMesh.nInteriorElements, Fe);


    return F;
}


VectorField cosYForce(double A, double k, const PolyMesh& theMesh) {

    // Forcing term vector initialization
    VectorField F;
    F.assign(theMesh.nInteriorElements, {0,0,0});


    // Forcing term computation
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        F[i].y = A * cos(k * theMesh.elements[i].centroid.x);
    }


    return F;
}