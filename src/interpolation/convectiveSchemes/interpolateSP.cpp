//
// Created by ruben on 4/12/23.
//

#include "convectiveSchemesOrthogonal.h"


GeometricVector interpolateSP(const GeometricVector& PhiC, const GeometricVector& PhiD) {

    return 0.5*(PhiC + PhiD);
}