//
// Created by ruben on 21/01/24.
//

#ifndef SOLVERNS_TENSOR_H
#define SOLVERNS_TENSOR_H

#include <vector>


class Tensor : public std::vector< std::vector<double> > {
public:
    // Tensor constructor and destructor
    Tensor();
    ~Tensor();
};


#endif //SOLVERNS_TENSOR_H
