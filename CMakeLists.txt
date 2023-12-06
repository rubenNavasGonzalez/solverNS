cmake_minimum_required(VERSION 3.22)

project(flowBetweenFlatPlates
        VERSION 0.0.1
        DESCRIPTION "Incompressible 3D laminar flow passed two parallel flat plates."
        LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

if (NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release)
endif ()

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/app)
#SET(CMAKE_CXX_FLAGS "-O3")

add_executable(flowBetweenFlatPlates app/main.cpp
        src/mesh/node/Node.cpp
        src/mesh/node/Node.h
        src/mesh/face/Face.cpp
        src/mesh/face/Face.h
        src/math/geometricVector/GeometricVector.cpp
        src/math/geometricVector/GeometricVector.h
        src/mesh/element/Element.cpp
        src/mesh/element/Element.h
        src/mesh/polyMesh/PolyMesh.cpp
        src/mesh/polyMesh/PolyMesh.h
        src/mesh/boundary/Boundary.cpp
        src/mesh/boundary/Boundary.h
        src/mesh/polyMesh/generatePolyMesh.cpp
        src/mesh/polyMesh/iterators2ElementVertices.cpp
        src/mesh/polyMesh/iterators2ElementVertices.h
        src/mesh/polyMesh/computeHexahedronVolume.cpp
        src/mesh/polyMesh/computeHexahedronVolume.h
        src/mesh/polyMesh/computeHexahedronCentroid.cpp
        src/mesh/polyMesh/computeHexahedronCentroid.h
        src/mesh/polyMesh/computeRectangleArea.cpp
        src/mesh/polyMesh/computeRectangleArea.h
        src/mesh/polyMesh/computeRectangleCentroid.cpp
        src/mesh/polyMesh/computeRectangleCentroid.h
        src/mesh/polyMesh/generateBoundaryMesh.cpp
        src/mesh/polyMesh/writeMesh2VTK.cpp
        src/fields/vectorField/VectorField.cpp
        src/fields/vectorField/VectorField.h
        src/fields/vectorField/initialize.cpp
        src/fields/scalarField/ScalarField.cpp
        src/fields/scalarField/ScalarField.h
        src/fields/scalarField/initialize.cpp
        src/boundaryConditions/vectorBoundaryConditions/VectorBoundaryConditions.cpp
        src/boundaryConditions/vectorBoundaryConditions/VectorBoundaryConditions.h
        src/boundaryConditions/vectorBoundaryConditions/addBC.cpp
        src/boundaryConditions/scalarBoundaryConditions/ScalarBoundaryConditions.cpp
        src/boundaryConditions/scalarBoundaryConditions/ScalarBoundaryConditions.h
        src/boundaryConditions/scalarBoundaryConditions/addBC.cpp
        src/finiteVolumeMethod/fvc/convectiveOrthogonal.cpp
        src/finiteVolumeMethod/fvc/convectiveOrthogonal.h
        src/finiteVolumeMethod/fvc/laplacianOrthogonal.cpp
        src/finiteVolumeMethod/fvc/laplacianOrthogonal.h
        src/interpolation/convectiveSchemes/convectiveSchemesOrthogonal.cpp
        src/interpolation/convectiveSchemes/convectiveSchemesOrthogonal.h
        src/interpolation/convectiveSchemes/interpolateUDS.cpp
        src/interpolation/convectiveSchemes/interpolateDDS.cpp
        src/fields/vectorField/applyBCs.cpp
        src/interpolation/convectiveSchemes/interpolateCDS_Orthogonal.cpp
        src/interpolation/convectiveSchemes/interpolateSUDS_Orthogonal.cpp
        src/interpolation/convectiveSchemes/interpolateQUICK_Orthogonal.cpp
        src/interpolation/convectiveSchemes/interpolateSMART_Orthogonal.cpp
        src/interpolation/convectiveSchemes/interpolateSP.cpp
        src/finiteVolumeMethod/fvc/divergence.cpp
        src/finiteVolumeMethod/fvc/divergence.h
        src/fields/scalarField/max.cpp
        src/math/sparseMatrix/SparseMatrix.cpp
        src/math/sparseMatrix/SparseMatrix.h
        src/math/sparseMatrix/getValue.cpp
        src/math/sparseMatrix/addValue.cpp
        src/finiteVolumeMethod/fvScalarEquation/FvScalarEquation.cpp
        src/finiteVolumeMethod/fvScalarEquation/FvScalarEquation.h
        src/finiteVolumeMethod/fvm/laplacianOrthogonal.cpp
        src/finiteVolumeMethod/fvm/laplacianOrthogonal.h
        src/finiteVolumeMethod/fvScalarEquation/constrain.cpp
        src/finiteVolumeMethod/fvScalarEquation/solve.cpp
        src/math/linearSolver/BiCGSTAB/solveBiCGSTAB.cpp
        src/math/linearSolver/BiCGSTAB/solveBiCGSTAB.h
        src/math/linearSolver/linearSolverConfig/LinearSolverConfig.cpp
        src/math/linearSolver/linearSolverConfig/LinearSolverConfig.h
)