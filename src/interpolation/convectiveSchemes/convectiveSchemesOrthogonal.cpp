//
// Created by ruben on 29/11/23.
//

#include "convectiveSchemesOrthogonal.h"


GeometricVector convectiveSchemesOrthogonal(const GeometricVector& PhiOwner, double xOwner, const GeometricVector& PhiOwnerFar, double xOwnerFar,
                                            const GeometricVector& PhiNeighbour, double xNeighbour, const GeometricVector& PhiNeighbourFar,
                                            double xNeighbourFar, double xF, double mDot, const std::string& scheme) {

    GeometricVector PhiF, PhiC, PhiU, PhiD;
    double xC, xU, xD;

    if (mDot >= 0) {

        xC = xOwner;
        xU = xOwnerFar;
        xD = xNeighbour;

        PhiC = PhiOwner;
        PhiU = PhiOwnerFar;
        PhiD = PhiNeighbour;
    } else {

        xC = xNeighbour;
        xU = xNeighbourFar;
        xD = xOwner;

        PhiC = PhiNeighbour;
        PhiU = PhiNeighbourFar;
        PhiD = PhiOwner;
    }


    if (scheme == "UDS") {

        PhiF = interpolateUDS(PhiC);
    } else if (scheme == "DDS") {

        PhiF = interpolateDDS(PhiD);
    } else if (scheme == "CDS") {

        PhiF = PhiC + (xF - xC)/(xD - xC)*(PhiD - PhiC);
    } else if (scheme == "SUDS") {

        PhiF = PhiU + (xF - xU)/(xC - xU)*(PhiC - PhiU);
    } else if (scheme == "QUICK") {

    } else if (scheme == "SMART") {

    } else {

        printf("ERROR. Invalid convective scheme selected!!\n");
    }

    return PhiF;
}