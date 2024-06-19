//
// Created by ruben on 26/12/23.
//

#include "fvc.h"


VectorField fvc::forcingTerm(const GeometricVector& Fe, const PolyMesh& theMesh) {

    // Forcing term vector initialization
    VectorField F;
    F.assign(theMesh.nInteriorElements, Fe);


    return F;
}