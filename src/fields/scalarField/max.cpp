//
// Created by ruben on 4/12/23.
//

#include <algorithm>
#include "ScalarField.h"


double ScalarField::max() {

    return *std::max_element(begin(), end());
}