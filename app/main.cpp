#include <cmath>
#include "../src/mesh/polyMesh/PolyMesh.h"
#include "../src/fields/vectorField/VectorField.h"
#include "../src/fields/scalarField/ScalarField.h"
#include "../src/boundaryConditions/vectorBoundaryConditions/VectorBoundaryConditions.h"
#include "../src/boundaryConditions/scalarBoundaryConditions/ScalarBoundaryConditions.h"
#include "../src/finiteVolumeMethod/fvc/laplacian.h"

//Mesh parameters
double Lx = 1, Ly = 1, Lz = 1;
int Nx = 16, Ny = 16, Nz = 1;
double sx = 0, sy = 0, sz = 0;


int main() {

    // Generation of the mesh
    PolyMesh theMesh;
    theMesh.generatePolyMesh(Lx, Ly, Lz, Nx, Ny, Nz, sx, sy, sz);
    theMesh.generateBoundaryMesh(Nx, Ny, Nz);
    //theMesh.writeMesh2VTK("simResults");


    /*// Preallocation of the fields
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


    // Set boundary conditions
    VectorBoundaryConditions uBCs;
    uBCs.addBC("fixedValue", {1,0,0});
    uBCs.addBC("zeroGradient", {0,0,0});
    uBCs.addBC("fixedValue", {0,0,0});
    uBCs.addBC("fixedValue", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});

    ScalarBoundaryConditions pBCs;
    pBCs.addBC("zeroGradient", 0);
    pBCs.addBC("fixedValue", 0);
    pBCs.addBC("zeroGradient", 0);
    pBCs.addBC("zeroGradient", 0);
    pBCs.addBC("periodic", 0);
    pBCs.addBC("periodic", 0);*/


    // Verification of the diffusive term
    double laplacianUDifferenceX, laplacianUDifferenceY;
    double maxX{}, maxY{};

    VectorField u;
    u.initialize(theMesh.nElements);

    VectorBoundaryConditions uBCs;
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("empty", {0,0,0});
    uBCs.addBC("empty", {0,0,0});

    for (int i = 0; i < theMesh.nElements; ++i) {

        u.field[i].x = cos(2*M_PI*theMesh.elements[i].centroid.x)*sin(2*M_PI*theMesh.elements[i].centroid.y);
        u.field[i].y = -sin(2*M_PI*theMesh.elements[i].centroid.x)*cos(2*M_PI*theMesh.elements[i].centroid.y);
        u.field[i].z = 0;
    }

    VectorField laplacianU = fvc::laplacian(u, theMesh, uBCs);

    for (int i = 0; i < theMesh.nElements; ++i) {

        laplacianUDifferenceX = fabs(laplacianU.field[i].x - (-8*pow(M_PI,2)*cos(2*M_PI*theMesh.elements[i].centroid.x)*sin(2*M_PI*theMesh.elements[i].centroid.y)));
        laplacianUDifferenceY = fabs(laplacianU.field[i].y - 8*pow(M_PI,2)*sin(2*M_PI*theMesh.elements[i].centroid.x)*cos(2*M_PI*theMesh.elements[i].centroid.y));

        if (laplacianUDifferenceX > maxX) {
            maxX = laplacianUDifferenceX;
        }
        if (laplacianUDifferenceY > maxY) {
            maxY = laplacianUDifferenceY;
        }
    }

    printf("The element size is: %f\n", Lx/Nx);
    printf("The maximum x error is: %f\n", maxX);
    printf("The maximum y error is: %f\n", maxY);

    return 0;
}
