//
// Created by ruben on 21/11/23.
//

#include "computeHexahedronVolume.h"


double computeHexahedronVolume(double x1, double x2, double y1, double y2, double z1, double z2) {

    return (x2 - x1)*(y2 - y1)*(z2 - z1);
}
