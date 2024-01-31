//
// Created by ruben on 31/01/24.
//

#include <cmath>
#include "computeQCrit.h"


ScalarField computeQCrit(const PolyMesh& theMesh, TensorField Phi) {

    // Auxiliary variables declaration
    ScalarField QCrit;
    QCrit.assign(theMesh.nInteriorElements, 0);


    // Loop for all the interior elements
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        QCrit[i] = 0.5*( pow(Phi[i].skewSymmetric().frobeniusNorm() ,2) - pow(Phi[i].symmetric().frobeniusNorm() ,2) );
    }


    return QCrit;
}