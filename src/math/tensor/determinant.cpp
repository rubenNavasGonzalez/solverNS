//
// Created by ruben on 8/02/24.
//

#include "Tensor.h"


double Tensor::determiant() {

    // Positive contribution computation
    double detPlus = this->at(0).at(0) * this->at(1).at(1) * this->at(2).at(2)
                     + this->at(0).at(1) * this->at(1).at(2) * this->at(2).at(0)
                     + this->at(1).at(0) * this->at(2).at(1) * this->at(0).at(2);


    // Negative contribution computation
    double detMinus = this->at(0).at(2) * this->at(1).at(1) * this->at(2).at(0)
                      + this->at(1).at(2) * this->at(2).at(1) * this->at(0).at(0)
                      + this->at(1).at(0) * this->at(0).at(1) * this->at(2).at(2);


    return detPlus - detMinus;
}