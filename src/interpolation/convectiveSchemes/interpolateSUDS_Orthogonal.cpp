//
// Created by ruben on 4/12/23.
//

#include "convectiveSchemesOrthogonal.h"


GeometricVector interpolateSUDS_Orthogonal(const GeometricVector& PhiU, const GeometricVector& PhiC, double xF, double xU, double xC) {

    return PhiU + (xF - xU)/(xC - xU)*(PhiC - PhiU);
}