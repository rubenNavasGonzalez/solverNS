//
// Created by ruben on 6/12/23.
//

#include "LinearSolverConfig.h"


// LinearSolverConfig default constructor and destructor
LinearSolverConfig::LinearSolverConfig() = default;
LinearSolverConfig::~LinearSolverConfig() = default;


// LinearSolverConfig constructor
LinearSolverConfig::LinearSolverConfig(std::string _solver, double _tolerance, int _maxIter) {

    solver = _solver;
    tolerance = _tolerance;
    maxIter = _maxIter;
}
