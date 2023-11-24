//
// Created by ruben on 19/11/23.
//

#include "Node.h"


// Node default constructor and destructor
Node::Node() = default;
Node::~Node() = default;

// Node constructor from a coordinate
Node::Node(double x_, double y_, double z_) {
    x = x_;
    y = y_;
    z = z_;
}


// Sum of nodes method
Node operator+(const Node &p1, const Node &p2) {
    return { p1.x + p2.x, p1.y + p2.y, p1.z + p2.z };
}

// Division of a node by a double method
Node operator/(const Node &p, const double &k) {
    return { p.x / k, p.y / k, p.z / k };
}

// Subtraction of nodes method
GeometricVector operator-(const Node &p1, const Node &p2) {
    return { p1.x - p2.x, p1.y - p2.y, p1.z - p2.z };
}
