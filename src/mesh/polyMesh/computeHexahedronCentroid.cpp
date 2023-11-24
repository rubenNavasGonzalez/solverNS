//
// Created by ruben on 21/11/23.
//

#include "computeHexahedronCentroid.h"


Node computeHexahedronCentroid(Node p1, Node p2, Node p3, Node p4, Node p5, Node p6, Node p7, Node p8) {

    return (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8)/8;
}