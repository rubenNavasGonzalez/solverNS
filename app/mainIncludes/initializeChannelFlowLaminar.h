//
// Created by ruben on 29/12/23.
//

#ifndef SOLVERNS_INITIALIZECHANNELFLOWLAMINAR_H
#define SOLVERNS_INITIALIZECHANNELFLOWLAMINAR_H


for (int i = 0; i < theMesh.nInteriorElements; ++i) {


    if (theMesh.elements[i].centroid.x == theMesh.elements[0].centroid.x) {

        u[i] = {1,0,0};
    } else {

        u[i] = {0,0,0};
    }
}

#endif //SOLVERNS_INITIALIZECHANNELFLOWLAMINAR_H
