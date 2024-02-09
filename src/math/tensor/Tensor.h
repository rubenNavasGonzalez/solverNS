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
    double trace();
    double determiant();
    double invariantP();
    double invariantQ();
    double invariantR();
    double invariantV2();
    friend Tensor operator+(const Tensor& T1, const Tensor& T2);
    friend Tensor operator-(const Tensor& T1, const Tensor& T2);
    friend Tensor operator*(const double& k, const Tensor& T);
    friend Tensor operator*(const Tensor& T1, const Tensor& T2);
};

#endif //SOLVERNS_TENSOR_H