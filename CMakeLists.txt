cmake_minimum_required(VERSION 3.22)

project(solverNS
        VERSION 0.0.1
        DESCRIPTION "Incompressible Navier Stokes solver using the FSM algorithm."
        LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

if (NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release)
endif ()

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/app)
#SET(CMAKE_CXX_FLAGS "-O3")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

#find_package(OpenMP REQUIRED)
#SET(CMAKE_CXX_FLAGS "-fopenmp")

add_executable(solverNS app/main.cpp
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
        src/fields/scalarField/ScalarField.cpp
        src/fields/scalarField/ScalarField.h
        src/boundaryConditions/vectorBoundaryConditions/VectorBoundaryConditions.cpp
        src/boundaryConditions/vectorBoundaryConditions/VectorBoundaryConditions.h
        src/boundaryConditions/vectorBoundaryConditions/addBC.cpp
        src/boundaryConditions/scalarBoundaryConditions/ScalarBoundaryConditions.cpp
        src/boundaryConditions/scalarBoundaryConditions/ScalarBoundaryConditions.h
        src/boundaryConditions/scalarBoundaryConditions/addBC.cpp
        src/finiteVolumeMethod/fvc/convective.cpp
        src/finiteVolumeMethod/fvc/laplacianOrthogonal.cpp
        src/finiteVolumeMethod/fvc/divergence.cpp
        src/fields/scalarField/max.cpp
        src/math/sparseMatrix/SparseMatrix.cpp
        src/math/sparseMatrix/SparseMatrix.h
        src/math/sparseMatrix/getValue.cpp
        src/math/sparseMatrix/addValue.cpp
        src/finiteVolumeMethod/fvScalarEquation/FvScalarEquation.cpp
        src/finiteVolumeMethod/fvScalarEquation/FvScalarEquation.h
        src/finiteVolumeMethod/fvm/laplacianOrthogonal.cpp
        src/finiteVolumeMethod/fvScalarEquation/constrain.cpp
        src/finiteVolumeMethod/fvScalarEquation/solve.cpp
        src/math/linearSolver/BiCGSTAB/solveBiCGSTAB.cpp
        src/math/linearSolver/BiCGSTAB/solveBiCGSTAB.h
        src/math/linearSolver/linearSolverConfig/LinearSolverConfig.cpp
        src/math/linearSolver/linearSolverConfig/LinearSolverConfig.h
        src/fields/scalarField/maxAbs.cpp
        src/finiteVolumeMethod/fvc/gradient.cpp
        src/finiteVolumeMethod/fvm/fvm.h
        src/finiteVolumeMethod/fvc/fvc.h
        src/interpolation/interpolateMDotFromElements2Faces/RhieChowInterpolation.cpp
        src/interpolation/interpolateMDotFromElements2Faces/RhieChowInterpolation.h
        src/math/linearSolver/SOR/solveSOR.cpp
        src/math/linearSolver/SOR/solveSOR.h
        src/interpolation/temporalAdvancement/computeTimeStepOrthogonal.cpp
        src/interpolation/temporalAdvancement/computeTimeStepOrthogonal.h
        src/fields/scalarField/min.cpp
        src/interpolation/interpolateFromElements2Nodes/interpolateVectorFieldFromElements2Nodes.cpp
        src/interpolation/interpolateFromElements2Nodes/interpolateVectorFieldFromElements2Nodes.h
        src/fields/vectorField/writeVectorField2VTK.cpp
        src/interpolation/interpolateFromElements2Nodes/interpolateScalarFieldFromElements2Nodes.cpp
        src/interpolation/interpolateFromElements2Nodes/interpolateScalarFieldFromElements2Nodes.h
        src/fields/scalarField/writeScalarField2VTK.cpp
        src/fields/vectorField/maxAbs.cpp
        src/finiteVolumeMethod/fvc/forcingTerm.cpp
        src/functionObjects/computeBulkVelocity/computeBulkVelocity.cpp
        src/functionObjects/computeBulkVelocity/computeBulkVelocity.h
        src/math/tensor/Tensor.cpp
        src/math/tensor/Tensor.h
        src/fields/tensorField/TensorField.cpp
        src/fields/tensorField/TensorField.h
        src/finiteVolumeMethod/fvc/curl.cpp
        src/IO/out/turbulentChannelFlow/writeTurbulentChannelFlowData2CSV.cpp
        src/IO/out/turbulentChannelFlow/writeTurbulentChannelFlowData2CSV.h
        src/IO/out/turbulentChannelFlow/writeTurbulentChannelFlowData2VTK.cpp
        src/IO/out/turbulentChannelFlow/writeTurbulentChannelFlowData2VTK.h
        src/math/linearSolver/CGS/solveCGS.cpp
        src/math/linearSolver/CGS/solveCGS.h
        src/finiteVolumeMethod/fvc/divergenceNVS.cpp
        src/fields/scalarField/rms.cpp
        src/math/tensor/transpose.cpp
        src/math/tensor/symmetric.cpp
        src/math/tensor/skewSymmetric.cpp
        src/math/tensor/frobeniusNorm.cpp
        src/functionObjects/computeQCrit/computeQCrit.cpp
        src/functionObjects/computeQCrit/computeQCrit.h
        src/IO/out/turbulentChannelFlow/writeTurbulentChannelFlowSetUp.cpp
        src/IO/out/turbulentChannelFlow/writeTurbulentChannelFlowSetUp.h
        src/math/linearSolver/CG/solveCG.cpp
        src/math/linearSolver/CG/solveCG.h
        src/fields/vectorField/initialize.cpp
        src/fields/vectorField/initializeTurbulentChannelFlowReTau180.cpp
        src/fields/vectorField/initializeTurbulentChannelFlowReTau180.h
        src/IO/in/readCSV.cpp
        src/IO/in/readCSV.h
        src/fields/vectorField/initializeFromSelectedTime.cpp
        src/fields/vectorField/initializeFromSelectedTime.h
        src/finiteVolumeMethod/fvScalarEquation/perturb.cpp
        src/functionObjects/probe/Probe.h
        src/functionObjects/probe/Probe.cpp
        src/functionObjects/probe/assign.cpp
        src/functionObjects/probe/writeField.cpp
        src/fields/vectorField/interpolateFromDifferentMesh.cpp
        src/fields/vectorField/interpolateFromDifferentMesh.h
        src/turbulenceModeling/computeTurbulentViscosity/computeTurbulentViscosity.cpp
        src/turbulenceModeling/computeTurbulentViscosity/computeTurbulentViscosity.h
        src/math/tensor/trace.cpp
        src/math/tensor/determinant.cpp
        src/turbulenceModeling/LES/modelSmagorinsky/modelSmagorinsky.cpp
        src/turbulenceModeling/LES/modelSmagorinsky/modelSmagorinsky.h
        src/math/tensor/invariantQ.cpp
        src/math/tensor/invariantR.cpp
        src/math/tensor/invariantV2.cpp
        src/turbulenceModeling/LES/modelWALE/modelWALE.cpp
        src/turbulenceModeling/LES/modelWALE/modelWALE.h
        src/turbulenceModeling/LES/modelVreman/modelVreman.cpp
        src/turbulenceModeling/LES/modelVreman/modelVreman.h
        src/turbulenceModeling/LES/modelVerstappen/modelVerstappen.cpp
        src/turbulenceModeling/LES/modelVerstappen/modelVerstappen.h
        src/turbulenceModeling/LES/modelLES.h
        src/turbulenceModeling/LES/modelS3PQR/modelS3PQR.cpp
        src/turbulenceModeling/LES/modelS3PQR/modelS3PQR.h
        src/math/tensor/invariantP.cpp
        src/finiteVolumeMethod/fvScalarEquation/changeSign.cpp
)

#target_compile_options(solverNS PRIVATE ${OpenMP_CXX_FLAGS})
#target_link_libraries(solverNS PRIVATE OpenMP::OpenMP_CXX)
