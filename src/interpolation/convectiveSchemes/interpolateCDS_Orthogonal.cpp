//
// Created by ruben on 4/12/23.
//

#include "convectiveSchemesOrthogonal.h"


GeometricVector interpolateCDS_Orthogonal(const GeometricVector& PhiC, const GeometricVector& PhiD, double xF, double xC, double xD) {

    return PhiC + (xF - xC)/(xD - xC)*(PhiD - PhiC);
}