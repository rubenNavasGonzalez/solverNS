//
// Created by ruben on 19/11/23.
//

#ifndef FLOWBETWEENFLATPLATES_NODE_H
#define FLOWBETWEENFLATPLATES_NODE_H

#include "../../math/geometricVector/GeometricVector.h"


class Node {
public:
    // Node x,y,z coordinates
    double x{}, y{}, z{};


    // Node constructors and destructor
    Node();
    ~Node();
    Node(double x_, double y_, double z_);


    // Operator overloading methods
    friend Node operator+(const Node& p1, const Node& p2);
    friend Node operator/(const Node& p, const double& k);
    friend GeometricVector operator-(const Node& p1, const Node& p2);
};


#endif //FLOWBETWEENFLATPLATES_NODE_H
