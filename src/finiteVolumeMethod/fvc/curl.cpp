//
// Created by ruben on 21/01/24.
//

#include "fvc.h"


VectorField fvc::curl(const TensorField& gradPhi, const PolyMesh& theMesh) {

    // Initialize the curl field
    VectorField curl;
    curl.resize(theMesh.nInteriorElements);


    // Loop over all the interior elements and assemble curl vector field
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        curl[i].x = gradPhi[i][1][2] - gradPhi[i][2][1];
        curl[i].y = gradPhi[i][2][0] - gradPhi[i][0][2];
        curl[i].z = gradPhi[i][0][1] - gradPhi[i][1][0];
    }


    return curl;
}