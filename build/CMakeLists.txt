cmake_minimum_required(VERSION 3.22)
project(flowBetweenFlatPlates)

set(CMAKE_CXX_STANDARD 17)

add_executable(flowBetweenFlatPlates ../app/main.cpp
        ../src/mesh/node/Node.cpp
        ../src/mesh/node/Node.h
        ../src/mesh/face/Face.cpp
        ../src/mesh/face/Face.h
        ../src/math/geometricVector/GeometricVector.cpp
        ../src/math/geometricVector/GeometricVector.h
        ../src/mesh/element/Element.cpp
        ../src/mesh/element/Element.h
        ../src/mesh/polyMesh/PolyMesh.cpp
        ../src/mesh/polyMesh/PolyMesh.h
        ../src/mesh/boundary/Boundary.cpp
        ../src/mesh/boundary/Boundary.h
        ../src/mesh/polyMesh/generatePolyMesh.cpp
        ../src/mesh/polyMesh/iterators2ElementVertices.cpp
        ../src/mesh/polyMesh/iterators2ElementVertices.h
        ../src/mesh/polyMesh/computeHexahedronVolume.cpp
        ../src/mesh/polyMesh/computeHexahedronVolume.h
        ../src/mesh/polyMesh/computeHexahedronCentroid.cpp
        ../src/mesh/polyMesh/computeHexahedronCentroid.h
        ../src/mesh/polyMesh/computeRectangleArea.cpp
        ../src/mesh/polyMesh/computeRectangleArea.h
        ../src/mesh/polyMesh/computeRectangleCentroid.cpp
        ../src/mesh/polyMesh/computeRectangleCentroid.h
        ../src/mesh/polyMesh/generateBoundaryMesh.cpp
        ../src/mesh/polyMesh/writeMesh2VTK.cpp)
