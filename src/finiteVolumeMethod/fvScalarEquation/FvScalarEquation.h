//
// Created by ruben on 6/12/23.
//

#ifndef FLOWBETWEENFLATPLATES_FVSCALAREQUATION_H
#define FLOWBETWEENFLATPLATES_FVSCALAREQUATION_H

#include "../../math/sparseMatrix/SparseMatrix.h"
#include "../../fields/scalarField/ScalarField.h"
#include "../../mesh/polyMesh/PolyMesh.h"
#include "../../boundaryConditions/scalarBoundaryConditions/ScalarBoundaryConditions.h"
#include "../../math/linearSolver/linearSolverConfig/LinearSolverConfig.h"


class FvScalarEquation {
public:
    // FvScalarEquation attributes
    SparseMatrix A;
    ScalarField b;


    // FvScalarEquation constructor and destructor
    FvScalarEquation();
    ~FvScalarEquation();


    // FvScalarEquation methods
    void constrain(const PolyMesh& theMesh, const ScalarBoundaryConditions& PhiBCs);
    ScalarField solve(const LinearSolverConfig& theLinearSolverConfig, const ScalarField& PhiOld);
    void perturb();
    void changeSign();
};


// Overload operator ==
FvScalarEquation operator==(const SparseMatrix& A, const ScalarField& b);

#endif //FLOWBETWEENFLATPLATES_FVSCALAREQUATION_H
