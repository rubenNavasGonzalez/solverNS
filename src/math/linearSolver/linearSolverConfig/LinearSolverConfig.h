//
// Created by ruben on 6/12/23.
//

#ifndef FLOWBETWEENFLATPLATES_LINEARSOLVERCONFIG_H
#define FLOWBETWEENFLATPLATES_LINEARSOLVERCONFIG_H

#include <string>


class LinearSolverConfig {
public:
    // Linear solver attributes
    std::string solver;
    std::string residualNorm;
    double tolerance{};
    int maxIter{};


    // LinearSolverConfig constructors and destructor
    LinearSolverConfig();
    LinearSolverConfig(std::string _solver, std::string _residualNorm, double _tolerance, int _maxIter);
    ~LinearSolverConfig();
};


#endif //FLOWBETWEENFLATPLATES_LINEARSOLVERCONFIG_H
