//
// Created by ruben on 15/12/23.
//

#include <algorithm>
#include "ScalarField.h"


double ScalarField::min() {

    return *std::min_element(begin(), end());
}