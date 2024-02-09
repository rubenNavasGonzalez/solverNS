//
// Created by ruben on 7/02/24.
//

#include "Probe.h"


// Probe default constructor and destructor
Probe::Probe() = default;
Probe::~Probe() = default;


// Probe overloaded constructor
Probe::Probe(double x, double y, double z) {

    this->p = {x, y, z};

    this->iElements.assign(2,0);
    this->d.assign(2,0);
}
