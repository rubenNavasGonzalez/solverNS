//
// Created by ruben on 10/04/24.
//

#include "ScalarField.h"
#include "initializeFromSelectedTime.h"
#include "initializeFromDifferentMesh.h"


void ScalarField::initialize(const PolyMesh& theMesh, double t, int mode) {

    if (t == 0) {

        // Initialize the velocity field according to the selected case
        this->assign(theMesh.nInteriorElements, 0);

    } else {

        if (mode == 0) {

            // Initialize the velocity field from time t and equal mesh
            *this = initializeScalarFieldFromSelectedTime(theMesh, t);

        } else if (mode == 1) {

            // Initialize the velocity field from time t and different mesh
            *this = interpolateScalarFieldFromDifferentMesh(theMesh, t);

        } else {

            printf("Error. No correct initialization mode selected. Valid ones are 0 (start from time t) and 1 "
                   "(interpolate from different mesh at time t). \n");
            std::exit(EXIT_FAILURE);
        }
    }
}