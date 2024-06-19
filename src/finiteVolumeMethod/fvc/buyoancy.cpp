//
// Created by ruben on 5/04/24.
//

#include "fvc.h"


VectorField fvc::buyoancy(const PolyMesh& theMesh, double Pr, const ScalarField &T, const GeometricVector& e) {

    // Preallocate the buoyancy field
    VectorField buyoancy;
    buyoancy.assign(theMesh.nInteriorElements, {0,0,0});


    // Loop over all the interior elements
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        buyoancy[i].x = Pr * T[i] * e.x;
        buyoancy[i].y = Pr * T[i] * e.y;
        buyoancy[i].z = Pr * T[i] * e.z;
    }


    return buyoancy;
}