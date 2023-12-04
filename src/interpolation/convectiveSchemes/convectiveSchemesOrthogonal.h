//
// Created by ruben on 29/11/23.
//

#ifndef FLOWBETWEENFLATPLATES_CONVECTIVESCHEMESORTHOGONAL_H
#define FLOWBETWEENFLATPLATES_CONVECTIVESCHEMESORTHOGONAL_H

#include <string>
#include "../../math/geometricVector/GeometricVector.h"
#include "../../mesh/node/Node.h"


GeometricVector convectiveSchemesOrthogonal(const GeometricVector& PhiOwner, const Node& xOwner, const GeometricVector& PhiOwnerFar, const Node& xOwnerFar,
                                            const GeometricVector& PhiNeighbour, const Node& xNeighbour, const GeometricVector& PhiNeighbourFar,
                                            const Node& xNeighbourFar, const Node& xF, double mDot, const GeometricVector& Sf, const std::string& scheme);

GeometricVector interpolateUDS(const GeometricVector& PhiU);
GeometricVector interpolateDDS(const GeometricVector& PhiD);
GeometricVector interpolateCDS_Orthogonal(const GeometricVector& PhiC, const GeometricVector& PhiD, double xF, double xC, double xD);
GeometricVector interpolateSUDS_Orthogonal(const GeometricVector& PhiU, const GeometricVector& PhiC, double xF, double xU, double xC);
GeometricVector interpolateQUICK_Orthogonal(const GeometricVector& PhiU, const GeometricVector& PhiC, const GeometricVector& PhiD, double xF, double xU, double xC, double xD);
GeometricVector interpolateSMART_Orthogonal(const GeometricVector& PhiU, const GeometricVector& PhiC, const GeometricVector& PhiD, double xF, double xU, double xC, double xD);
GeometricVector interpolateSP(const GeometricVector& PhiC, const GeometricVector& PhiD);


#endif //FLOWBETWEENFLATPLATES_CONVECTIVESCHEMESORTHOGONAL_H
