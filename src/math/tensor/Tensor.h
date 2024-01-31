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


    // Tensor methods
    Tensor transpose();
    Tensor symmetric();
    Tensor skewSymmetric();
    double frobeniusNorm();
    friend Tensor operator+(const Tensor& T1, const Tensor& T2);
    friend Tensor operator-(const Tensor& T1, const Tensor& T2);
    friend Tensor operator*(const double& k, const Tensor& T);
};

#endif //SOLVERNS_TENSOR_H