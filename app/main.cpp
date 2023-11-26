#include "../src/mesh/polyMesh/PolyMesh.h"
#include "../src/fields/vectorField/VectorField.h"
#include "../src/fields/scalarField/ScalarField.h"

//Mesh parameters
double Lx = 1, Ly = 1.5, Lz = 1;
int Nx = 3, Ny = 3, Nz = 2;
double sx = 0, sy = 0, sz = 0;


int main() {

    // Generation of the mesh
    PolyMesh theMesh;
    theMesh.generatePolyMesh(Lx, Ly, Lz, Nx, Ny, Nz, sx, sy, sz);
    theMesh.generateBoundaryMesh(Nx, Ny, Nz);
    //theMesh.writeMesh2VTK("simResults");


    // Preallocation of the fields
    VectorField u, uPred, uFuture;
    u.initialize(theMesh.nElements);
    uPred.initialize(theMesh.nElements);
    uFuture.initialize(theMesh.nElements);

    VectorField R, RPrev;
    R.initialize(theMesh.nElements);
    RPrev.initialize(theMesh.nElements);

    ScalarField mDot;
    mDot.initialize(theMesh.nFaces);

    ScalarField p;
    p.initialize(theMesh.nElements);


    return 0;
}
