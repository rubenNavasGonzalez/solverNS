#include <cmath>
#include <fstream>
#include "../src/interpolation/temporalAdvancement/computeTimeStepOrthogonal.h"
#include "../src/finiteVolumeMethod/fvc/fvc.h"
#include "../src/finiteVolumeMethod/fvm/fvm.h"
#include "../src/finiteVolumeMethod/fvScalarEquation/FvScalarEquation.h"
#include "../src/interpolation/interpolateMDotFromElements2Faces/RhieChowInterpolation.h"


//Mesh parameters
double Lx = 10, Ly = 1, Lz = 1;
int Nx = 32, Ny = 32, Nz = 1;
double sx = 0, sy = 0, sz = 0;


int main() {

    // Generation of the mesh
    PolyMesh theMesh;
    theMesh.generatePolyMesh(Lx, Ly, Lz, Nx, Ny, Nz, sx, sy, sz);
    theMesh.generateBoundaryMesh(Nx, Ny, Nz);
    theMesh.writeMesh2VTK("simResults");


    double t = 0;
    double DeltaT;
    double nu = 0.01;

    VectorField u;
    u.initialize(theMesh.nInteriorElements);

    VectorBoundaryConditions uBCs;
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("empty", {0,0,0});
    uBCs.addBC("empty", {0,0,0});

    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        u.field[i].x = cos(2*M_PI*theMesh.elements[i].centroid.x)*sin(2*M_PI*theMesh.elements[i].centroid.y);
        u.field[i].y = -sin(2*M_PI*theMesh.elements[i].centroid.x)*cos(2*M_PI*theMesh.elements[i].centroid.y);
        u.field[i].z = 0;
    }

    ScalarField p;
    p.initialize(theMesh.nInteriorElements);

    ScalarBoundaryConditions pBCs;
    pBCs.addBC("periodic", 0);
    pBCs.addBC("periodic", 0);
    pBCs.addBC("periodic", 0);
    pBCs.addBC("periodic", 0);
    pBCs.addBC("empty", 0);
    pBCs.addBC("empty", 0);

    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        p.field[i] = -0.5*(pow(cos(2*M_PI*theMesh.elements[i].centroid.x),2) + pow(cos(2*M_PI*theMesh.elements[i].centroid.y),2));
    }

    LinearSolverConfig pSolver("BiCGSTAB", "L2", 1e-9, 1e6);

    VectorField RPrev;
    VectorField uNodeValue;
    ScalarField DeltaPValue;

    uNodeValue.field.push_back(u[50]);
    DeltaPValue.field.push_back(p[50] - p[0]);

    while (t < 3) {

        DeltaT = computeTimeStepOrthogonal(theMesh, u, nu);
        t += DeltaT;
        printf("\nTime = %f s \n", t);

        ScalarField mDot = RhieChowInterpolation(u, p, theMesh, DeltaT, uBCs, pBCs);
        VectorField convU = fvc::convectiveOrthogonal(mDot, u, theMesh, uBCs, "CDS");
        VectorField diffU = fvc::laplacianOrthogonal(u, theMesh, uBCs);

        VectorField R = nu*diffU - convU;

        if (t == DeltaT) {

            RPrev = R;
        }

        VectorField uPred = u + (DeltaT/1)*(1.5*R - 0.5*RPrev);

        ScalarField divUPred = fvc::divergence(uPred, theMesh, uBCs);
        SparseMatrix laplacianMatrixU = fvm::laplacianOrthogonal(theMesh);

        FvScalarEquation pEqn  =  laplacianMatrixU == (1/DeltaT)*divUPred;
        pEqn.constrain(theMesh, pBCs);
        printf("\tSolving Poisson Equation...\n");
        ScalarField pNew = pEqn.solve(pSolver, p);
        DeltaPValue.field.push_back(pNew[50] - pNew[0]);

        VectorField gradP = fvc::gradient(pNew, theMesh, pBCs);

        VectorField uNew = uPred - (DeltaT/1)*gradP;
        uNodeValue.field.push_back(uNew[50]);

        ScalarField divUNew = fvc::divergence(uNew, theMesh, uBCs);
        printf("\n\tThe maximum value of divUNew is: %E \n", divUNew.maxAbs());

        u = uNew;
        p = pNew;
        RPrev = R;
    }


    std::ofstream output;
    output.open ("data.txt");
    for (int i = 0; i < uNodeValue.field.size(); ++i) {

        output << uNodeValue[i].x << "\t" << uNodeValue[i].y << "\t" << DeltaPValue[i] << "\n";
    }
    output.close();

    return 0;
}