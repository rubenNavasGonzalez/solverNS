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


    return laplacian;
}