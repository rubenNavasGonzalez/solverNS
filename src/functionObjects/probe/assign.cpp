//
// Created by ruben on 7/02/24.
//

#include "Probe.h"


void Probe::assign(const PolyMesh& theMesh) {

    // Auxiliary variables initialization
    double dMin1 = 1e24;
    double dMin2 = 1e25;
    double dMag;
    GeometricVector _d;


    // Loop for all the interior elements
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {


        // Distance vector between the element i and the probe
        _d = theMesh.elements[i].centroid - this->p;
        dMag = _d.mag();


        // Check if the distance is minimum
        if ( dMag < dMin2 ) {

            if ( dMag <= dMin1 ) {

                dMin2 = dMin1;
                dMin1 = dMag;

                this->iElements[1] =  this->iElements[0];
                this->iElements[0] = i;

            } else {

                dMin2 = dMag;
                this->iElements[1] = i;
            }
        }
    }

    this->d[0] = dMin1;
    this->d[1] = dMin2;
}