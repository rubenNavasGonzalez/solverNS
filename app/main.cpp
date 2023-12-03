#include <cmath>
#include "../src/mesh/polyMesh/PolyMesh.h"
#include "../src/fields/vectorField/VectorField.h"
#include "../src/fields/scalarField/ScalarField.h"
#include "../src/boundaryConditions/vectorBoundaryConditions/VectorBoundaryConditions.h"
#include "../src/boundaryConditions/scalarBoundaryConditions/ScalarBoundaryConditions.h"
#include "../src/finiteVolumeMethod/fvc/laplacianOrthogonal.h"
#include "../src/finiteVolumeMethod/fvc/convectiveOrthogonal.h"


//Mesh parameters
double Lx = 1, Ly = 1, Lz = 1;
int Nx = 3, Ny = 3, Nz = 2;
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
    double differenceX, differenceY, differenceZ;
    double maxX{}, maxY{}, maxZ{};

    VectorField u;

    VectorBoundaryConditions uBCs;
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});

    u.initialize(theMesh.nElements);

    for (int i = 0; i < theMesh.nElements; ++i) {

        u.field[i].x = cos(2*M_PI*theMesh.elements[i].centroid.x)*sin(2*M_PI*theMesh.elements[i].centroid.y)*cos(2*M_PI*theMesh.elements[i].centroid.z);
        u.field[i].y = -sin(2*M_PI*theMesh.elements[i].centroid.x)*cos(2*M_PI*theMesh.elements[i].centroid.y)*cos(2*M_PI*theMesh.elements[i].centroid.z);
        u.field[i].z = 0;
    }
    //u.applyBCs(theMesh, uBCs);

    VectorField laplacianU = fvc::laplacianOrthogonal(u, theMesh, uBCs);

    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        differenceX = fabs(laplacianU.field[i].x - (-12*pow(M_PI,2)*u.field[i].x));
        differenceY = fabs(laplacianU.field[i].y + 12*pow(M_PI,2)*u.field[i].y);
        differenceZ = fabs(laplacianU.field[i].z - 0);

        if (differenceX > maxX) {
            maxX = differenceX;
        }
        if (differenceY > maxY) {
            maxY = differenceY;
        }
        if (differenceZ > maxZ) {
            maxZ = differenceZ;
        }
    }

    printf("\nThe number of elements is: %ix%ix%i \n\n", Nx, Ny, Nz);
    printf("Diffusive term verification: \n");
    printf("\tThe maximum x error is: %f \n", maxX);
    printf("\tThe maximum y error is: %f \n", maxY);
    printf("\tThe maximum z error is: %f \n\n", maxZ);


    // Verification of the convective term
    maxX = 0; maxY = 0; maxZ = 0;

     for (int i = 0; i < theMesh.nElements; ++i) {

         u.field[i].x = cos(2*M_PI*theMesh.elements[i].centroid.x)*sin(2*M_PI*theMesh.elements[i].centroid.y)*cos(2*M_PI*theMesh.elements[i].centroid.z);
         u.field[i].y = -sin(2*M_PI*theMesh.elements[i].centroid.x)*cos(2*M_PI*theMesh.elements[i].centroid.y)*cos(2*M_PI*theMesh.elements[i].centroid.z);
         u.field[i].z = 0;
     }

     ScalarField mDot;
     VectorField convU = fvc::convectiveOrthogonal(mDot, u, theMesh, uBCs, "UDS");

     for (int i = 0; i < theMesh.nInteriorElements; ++i) {

         differenceX = fabs(convU.field[i].x + M_PI*sin(4*M_PI*theMesh.elements[i].centroid.x)*pow(cos(2*M_PI*theMesh.elements[i].centroid.z),2));
         differenceY = fabs(convU.field[i].y + M_PI*sin(4*M_PI*theMesh.elements[i].centroid.y)*pow(cos(2*M_PI*theMesh.elements[i].centroid.z),2));
         differenceZ = fabs(convU.field[i].z - 0);

         if (differenceX > maxX) {
             maxX = differenceX;
         }
         if (differenceY > maxY) {
             maxY = differenceY;
         }
         if (differenceZ > maxZ) {
             maxZ = differenceZ;
         }
     }

     printf("Convective term verification: \n\n");
     printf("\tThe maximum x error is: %f \n", maxX);
     printf("\tThe maximum y error is: %f\n", maxY);
     printf("\tThe maximum z error is: %f\n", maxZ);


    return 0;
}
