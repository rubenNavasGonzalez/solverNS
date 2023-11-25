#include "../src/mesh/polyMesh/PolyMesh.h"

//Mesh parameters
double Lx = 1, Ly = 1.5, Lz = 1;
int Nx = 3, Ny = 3, Nz = 2;
double sx = 0, sy = 0, sz = 0;


int main() {

    // Generation of the mesh
    PolyMesh theMesh;
    theMesh.generatePolyMesh(Lx, Ly, Lz, Nx, Ny, Nz, sx, sy, sz);
    theMesh.generateBoundaryMesh(Nx, Ny, Nz);
    theMesh.writeMesh2VTK("simResults");

    return 0;
}
