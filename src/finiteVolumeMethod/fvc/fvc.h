//
// Created by ruben on 7/12/23.
//

#ifndef FLOWBETWEENFLATPLATES_FVC_H
#define FLOWBETWEENFLATPLATES_FVC_H

#include "../../fields/vectorField/VectorField.h"
#include "../../fields/scalarField/ScalarField.h"
#include "../../fields/tensorField/TensorField.h"
#include "../../boundaryConditions/scalarBoundaryConditions/ScalarBoundaryConditions.h"


namespace fvc {

    VectorField convective(const ScalarField& mDot, const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs);
    ScalarField convective(const ScalarField& mDot, const ScalarField& Phi, const PolyMesh& theMesh, const ScalarBoundaryConditions& PhiBCs);

    ScalarField divergence(const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs);
    ScalarField divergenceNVS(const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs);

    VectorField gradient(const ScalarField& Phi, const PolyMesh& theMesh, const ScalarBoundaryConditions& PhiBCs);
    TensorField gradient(const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs);

    VectorField laplacianOrthogonal(const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs);
    ScalarField laplacianOrthogonal(const ScalarField& Phi, const PolyMesh& theMesh, const ScalarBoundaryConditions& PhiBCs);
    VectorField laplacianOrthogonal(const ScalarField& Gamma, const VectorField& Phi, const PolyMesh& theMesh, const VectorBoundaryConditions& PhiBCs);

    VectorField forcingTerm(const GeometricVector& Fe, const PolyMesh& theMesh);
    VectorField buyoancy(const PolyMesh& theMesh, double Pr, const ScalarField& T, const GeometricVector& e);

    VectorField curl(const TensorField& gradPhi, const PolyMesh& theMesh);
}

#endif //FLOWBETWEENFLATPLATES_FVC_H
