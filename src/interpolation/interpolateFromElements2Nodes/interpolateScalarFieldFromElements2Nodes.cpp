//
// Created by ruben on 18/12/23.
//

#include <vector>
#include "interpolateScalarFieldFromElements2Nodes.h"


ScalarField interpolateScalarFieldFromElements2Nodes(const ScalarField& Phi, const PolyMesh& theMesh, const ScalarBoundaryConditions& PhiBCs) {

    // Auxiliary variables
    int iNode;
    std::string BCType;
    double BCValue;


    // Preallocate the node-centered vector field
    ScalarField PhiNode;
    PhiNode.assign(theMesh.nNodes, 0);


    // Node frequency vector (used to compute the average)
    std::vector<int> nodeFrequency(theMesh.nNodes, 0);


    // Loop over all the interior elements
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        // Loop over all the nodes (j) of an interior element (i)
        for (int j = 0; j < theMesh.elements[i].iNodes.size(); ++j) {

            iNode = theMesh.elements[i].iNodes[j];

            PhiNode[iNode] += Phi[i];
            nodeFrequency[iNode] += 1;
        }
    }


    // Average the nodal field with the frequency vector
    for (int i = 0; i < PhiNode.size(); ++i) {

        PhiNode[i] /= nodeFrequency[i];
    }
    nodeFrequency.assign(theMesh.nNodes, 0);


    // Loop over all the boundaries
    for (int k = 0; k < theMesh.nBoundaries; ++k) {

        BCType = PhiBCs.type[k];
        BCValue = PhiBCs.value[k];


        // Check if the boundary condition is fixedValue type
        if (BCType == "fixedValue") {

            // Loop over all the faces of the boundary
            for (int i = theMesh.boundaries[k].startFace; i < theMesh.boundaries[k].startFace + theMesh.boundaries[k].nBoundaryFaces; ++i) {

                // Loop over all the nodes (j) of the boundary face (i) of the boundary (k)
                for (int j = 0; j < theMesh.faces[i].iNodes.size(); ++j) {

                    // Set the value of the node equal as the value of the boundary condition
                    iNode = theMesh.faces[i].iNodes[j];
                    PhiNode[iNode] = BCValue;
                }
            }
        }
    }


    return PhiNode;
}