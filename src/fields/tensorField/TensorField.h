//
// Created by ruben on 21/01/24.
//

#ifndef SOLVERNS_TENSORFIELD_H
#define SOLVERNS_TENSORFIELD_H

#include "../../math/tensor/Tensor.h"


class TensorField : public std::vector<Tensor> {
public:
    // TensorField constructor and destructor
    TensorField();
    ~TensorField();
};


#endif //SOLVERNS_TENSORFIELD_H
