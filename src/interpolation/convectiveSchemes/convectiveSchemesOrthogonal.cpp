//
// Created by ruben on 29/11/23.
//

#include "convectiveSchemesOrthogonal.h"


GeometricVector convectiveSchemesOrthogonal(const GeometricVector& PhiOwner, const Node& pOwner, const GeometricVector& PhiOwnerFar, const Node& pOwnerFar,
                                            const GeometricVector& PhiNeighbour, const Node& pNeighbour, const GeometricVector& PhiNeighbourFar,
                                            const Node& pNeighbourFar, const Node& pF, double mDot, const GeometricVector& Sf, const std::string& scheme) {

    GeometricVector PhiF, PhiC, PhiU, PhiD;
    double pOwnerValue, pOwnerFarValue, pNeighbourValue, pNeighbourFarValue, pFValue, xC, xU, xD, xF;

    if (Sf.x != 0) {

        pOwnerValue = pOwner.x;
        pOwnerFarValue = pOwnerFar.x;
        pNeighbourValue = pNeighbour.x;
        pNeighbourFarValue = pNeighbourFar.x;
        pFValue = pF.x;

    } else if (Sf.y != 0) {

        pOwnerValue = pOwner.y;
        pOwnerFarValue = pOwnerFar.y;
        pNeighbourValue = pNeighbour.y;
        pNeighbourFarValue = pNeighbourFar.y;
        pFValue = pF.y;

    } else {

        pOwnerValue = pOwner.z;
        pOwnerFarValue = pOwnerFar.z;
        pNeighbourValue = pNeighbour.z;
        pNeighbourFarValue = pNeighbourFar.z;
        pFValue = pF.z;
    }


    if (mDot >= 0) {

        xC = pOwnerValue;
        xU = pOwnerFarValue;
        xD = pNeighbourValue;
        xF = pFValue;

        PhiC = PhiOwner;
        PhiU = PhiOwnerFar;
        PhiD = PhiNeighbour;
    } else {

        xC = pNeighbourValue;
        xU = pNeighbourFarValue;
        xD = pOwnerValue;
        xF = pFValue;

        PhiC = PhiNeighbour;
        PhiU = PhiNeighbourFar;
        PhiD = PhiOwner;
    }


    if (scheme == "UDS") {

        PhiF = interpolateUDS(PhiC);
    } else if (scheme == "DDS") {

        PhiF = interpolateDDS(PhiD);
    } else if (scheme == "CDS") {

        PhiF = interpolateCDS_Orthogonal(PhiC, PhiD, xF, xC, xD);
    } else if (scheme == "SUDS") {

        PhiF = interpolateSUDS_Orthogonal(PhiU, PhiC, xF, xU, xC);
    } else if (scheme == "QUICK") {

        PhiF = interpolateQUICK_Orthogonal(PhiU, PhiC, PhiD, xF, xU, xC, xD);
    } else if (scheme == "SMART") {

        PhiF = interpolateSMART_Orthogonal(PhiU, PhiC, PhiD, xF, xU, xC, xD);
    } else if (scheme == "SP") {

        PhiF = interpolateSP(PhiC, PhiD);
    } else {

        printf("ERROR. Invalid convective scheme selected!!\n");
    }


    return PhiF;
}