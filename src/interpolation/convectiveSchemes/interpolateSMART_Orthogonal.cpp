//
// Created by ruben on 4/12/23.
//

#include "convectiveSchemesOrthogonal.h"


GeometricVector interpolateSMART_Orthogonal(const GeometricVector& PhiU, const GeometricVector& PhiC, const GeometricVector& PhiD, double xF, double xU, double xC, double xD) {

    double xFHat = (xF - xU)/(xD - xU);
    double xCHat = (xC - xU)/(xD - xU);
    GeometricVector PhiCHat = (PhiC - PhiU)/(PhiD - PhiU);

    GeometricVector PhiFHat, PhiF;

    // x component of Phi
    if (PhiCHat.x > 0 && PhiCHat.x < xCHat/3) {

        PhiFHat.x = - xFHat*(1 - 3*xCHat + 2*xFHat)/(xCHat*(xCHat - 1))*PhiCHat.x;
    } else if (PhiCHat.x > xCHat/3 && PhiCHat.x < xCHat/xFHat*(1 + xFHat - xCHat)) {

        PhiFHat.x = xFHat*(xFHat - xCHat)/(1 - xCHat) + xFHat*(xFHat - 1)/(xCHat*(xCHat - 1))*PhiCHat.x;
    } else if (PhiCHat.x > xCHat/xFHat*(1 + xFHat - xCHat) && PhiCHat.x < 1) {

        PhiFHat.x = 1;
    } else {

        PhiFHat.x = PhiCHat.x;
    }

    // y component of Phi
    if (PhiCHat.y > 0 && PhiCHat.y < xCHat/3) {

        PhiFHat.y = - xFHat*(1 - 3*xCHat + 2*xFHat)/(xCHat*(xCHat - 1))*PhiCHat.y;
    } else if (PhiCHat.y > xCHat/3 && PhiCHat.y < xCHat/xFHat*(1 + xFHat - xCHat)) {

        PhiFHat.y = xFHat*(xFHat - xCHat)/(1 - xCHat) + xFHat*(xFHat - 1)/(xCHat*(xCHat - 1))*PhiCHat.y;
    } else if (PhiCHat.y > xCHat/xFHat*(1 + xFHat - xCHat) && PhiCHat.y < 1) {

        PhiFHat.y = 1;
    } else {

        PhiFHat.y = PhiCHat.y;
    }

    // z component of Phi
    if (PhiCHat.z > 0 && PhiCHat.z < xCHat/3) {

        PhiFHat.z = - xFHat*(1 - 3*xCHat + 2*xFHat)/(xCHat*(xCHat - 1))*PhiCHat.z;
    } else if (PhiCHat.z > xCHat/3 && PhiCHat.z < xCHat/xFHat*(1 + xFHat - xCHat)) {

        PhiFHat.z = xFHat*(xFHat - xCHat)/(1 - xCHat) + xFHat*(xFHat - 1)/(xCHat*(xCHat - 1))*PhiCHat.z;
    } else if (PhiCHat.z > xCHat/xFHat*(1 + xFHat - xCHat) && PhiCHat.z < 1) {

        PhiFHat.z = 1;
    } else {

        PhiFHat.z = PhiCHat.z;
    }

    PhiF.x = PhiFHat.x*(PhiD.x - PhiU.x) + PhiU.x;
    PhiF.y = PhiFHat.y*(PhiD.y - PhiU.y) + PhiU.y;
    PhiF.z = PhiFHat.z*(PhiD.z - PhiU.z) + PhiU.z;


    return PhiF;
}
