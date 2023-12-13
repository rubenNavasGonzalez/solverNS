//
// Created by ruben on 6/12/23.
//

#include "fvm.h"


SparseMatrix fvm::laplacianOrthogonal(const PolyMesh& theMesh) {

    // Auxiliary variables definition
    double dONMag, SfMag, coeffValue;
    int iOwner, iNeighbour;


    // Preallocate the laplacian field
    SparseMatrix laplacian;
    /*laplacian.nRows = theMesh.nElements;
    laplacian.nCols = theMesh.nElements;*/
    laplacian.nRows = theMesh.nInteriorElements;
    laplacian.nCols = theMesh.nInteriorElements;
    laplacian.diagIndex.resize(theMesh.nInteriorElements);
    laplacian.diagValue.resize(theMesh.nInteriorElements);


    // Loop over all the interior faces
    for (int i = 0; i < theMesh.nInteriorFaces; ++i) {

        dONMag = theMesh.faces[i].dONMag;
        SfMag = theMesh.faces[i].SfMag;

        iOwner = theMesh.faces[i].iOwner;
        iNeighbour = theMesh.faces[i].iNeighbour;

        coeffValue = SfMag/dONMag;

        laplacian.diagIndex[iOwner] = {iOwner, iOwner};
        laplacian.diagValue[iOwner] -= coeffValue;
        laplacian.addValue(coeffValue, {iOwner, iNeighbour});

        laplacian.diagIndex[iNeighbour] = {iNeighbour, iNeighbour};
        laplacian.diagValue[iNeighbour] -= coeffValue;
        laplacian.addValue(coeffValue, {iNeighbour, iOwner});
    }


    // Average the matrix by the element volume
    for (int i = 0; i < laplacian.upperIndex.size(); ++i) {

        laplacian.upperValue[i] /= theMesh.elements[ laplacian.upperIndex[i][0] ].Vf;
    }
    for (int i = 0; i < laplacian.diagIndex.size(); ++i) {

        laplacian.diagValue[i] /= theMesh.elements[ laplacian.diagIndex[i][0] ].Vf;
    }
    for (int i = 0; i < laplacian.lowerIndex.size(); ++i) {

        laplacian.lowerValue[i] /= theMesh.elements[ laplacian.lowerIndex[i][0] ].Vf;
    }


    return laplacian;
}