//
// Created by ruben on 7/02/24.
//

#include "Probe.h"


// Probe default constructor and destructor
Probe::Probe() = default;
Probe::~Probe() = default;


// Probe overloaded constructor
Probe::Probe(double _x, double _y, double _z) {

    this->p = {_x, _y, _z};

    this->iElements.assign(2,0);
    this->d.assign(2,0);
}
