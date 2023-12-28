//
// Created by ruben on 26/12/23.
//

#include "fvc.h"


VectorField fvc::forcingTerm(const GeometricVector& Fe, const PolyMesh& theMesh) {

    // Forcing term vector initialization
    VectorField F;
    F.assign(theMesh.nInteriorElements, Fe);


/*
    // Average by the volume
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        F[i] /= theMesh.elements[i].Vf;
    }
*/


    return F;
}