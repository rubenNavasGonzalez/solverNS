//
// Created by ruben on 29/11/23.
//

#ifndef FLOWBETWEENFLATPLATES_CONVECTIVESCHEMESORTHOGONAL_H
#define FLOWBETWEENFLATPLATES_CONVECTIVESCHEMESORTHOGONAL_H

#include <string>
#include "../../math/geometricVector/GeometricVector.h"


GeometricVector convectiveSchemesOrthogonal(const GeometricVector& PhiOwner, double xOwner, const GeometricVector& PhiOwnerFar, double xOwnerFar,
                                            const GeometricVector& PhiNeighbour, double xNeighbour, const GeometricVector& PhiNeighbourFar,
                                            double xNeighbourFar, double xF, double mDot, const std::string& scheme);

#endif //FLOWBETWEENFLATPLATES_CONVECTIVESCHEMESORTHOGONAL_H
