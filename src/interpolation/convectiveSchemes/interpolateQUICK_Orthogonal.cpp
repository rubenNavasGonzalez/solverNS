//
// Created by ruben on 4/12/23.
//

#include "convectiveSchemesOrthogonal.h"


GeometricVector interpolateQUICK_Orthogonal(const GeometricVector& PhiU, const GeometricVector& PhiC, const GeometricVector& PhiD, double xF, double xU, double xC, double xD) {

    double xFHat = (xF - xU)/(xD - xU);
    double xCHat = (xC - xU)/(xD - xU);
    GeometricVector PhiCHat = (PhiC - PhiU)/(PhiD - PhiU);

    GeometricVector PhiFHat, PhiF;
    PhiFHat.x = xFHat + xFHat*(xFHat - 1)/(xCHat*(xCHat - 1))*(PhiCHat.x - xCHat);
    PhiFHat.y = xFHat + xFHat*(xFHat - 1)/(xCHat*(xCHat - 1))*(PhiCHat.y - xCHat);
    PhiFHat.z = xFHat + xFHat*(xFHat - 1)/(xCHat*(xCHat - 1))*(PhiCHat.z - xCHat);

    PhiF.x = PhiFHat.x*(PhiD.x - PhiU.x) + PhiU.x;
    PhiF.y = PhiFHat.y*(PhiD.y - PhiU.y) + PhiU.y;
    PhiF.z = PhiFHat.z*(PhiD.z - PhiU.z) + PhiU.z;


    return PhiF;
}