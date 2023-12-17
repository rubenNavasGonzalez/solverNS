#include <fstream>
#include "../src/interpolation/temporalAdvancement/computeTimeStepOrthogonal.h"
#include "../src/finiteVolumeMethod/fvc/fvc.h"
#include "../src/finiteVolumeMethod/fvm/fvm.h"
#include "../src/finiteVolumeMethod/fvScalarEquation/FvScalarEquation.h"
#include "../src/interpolation/interpolateMDotFromElements2Faces/RhieChowInterpolation.h"


// Mesh parameters
double Lx = 10, Ly = 1, Lz = 1;
int Nx = 128, Ny = 32, Nz = 8;
double sx = 0, sy = 0, sz = 0;


int main() {

    // Generation of the mesh
    PolyMesh theMesh;
    theMesh.generatePolyMesh(Lx, Ly, Lz, Nx, Ny, Nz, sx, sy, sz);
    theMesh.generateBoundaryMesh(Nx, Ny, Nz);
    theMesh.writeMesh2VTK("simResults");


    double t = 0;
    double DeltaT;
    double nu = 0.05;

    VectorField u;
    u.assign(theMesh.nInteriorElements, {0,0,0});

    VectorBoundaryConditions uBCs;
    uBCs.addBC("fixedValue", {1,0,0});
    uBCs.addBC("zeroGradient", {0,0,0});
    uBCs.addBC("fixedValue", {0,0,0});
    uBCs.addBC("fixedValue", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});

    ScalarField p;
    p.assign(theMesh.nInteriorElements, 0);

    ScalarBoundaryConditions pBCs;
    pBCs.addBC("zeroGradient", 0);
    pBCs.addBC("fixedValue", 0);
    pBCs.addBC("zeroGradient", 0);
    pBCs.addBC("zeroGradient", 0);
    pBCs.addBC("periodic", 0);
    pBCs.addBC("periodic", 0);


    LinearSolverConfig pSolver("BiCGSTAB", "L2", 1e-6, 1e6);

    VectorField RPrev;


    while (t < 10) {

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

        VectorField gradP = fvc::gradient(pNew, theMesh, pBCs);

        VectorField uNew = uPred - (DeltaT/1)*gradP;

        ScalarField divUNew = fvc::divergence(uNew, theMesh, uBCs);
        printf("\n\tThe maximum value of divUNew is: %E \n", divUNew.max());


        u = uNew;
        p = pNew;
        RPrev = R;
    }


    int NxHalf = Nx/2;
    int NyHalf = Ny/2;
    int NzHalf = Nz/2;

    std::ofstream output;
    output.open ("dataU.txt");
    for (int i = 0; i < Ny; ++i) {

        output << u[Nx*Ny*NzHalf + Nx*i + NxHalf].x << "\n";
    }
    output.close();

    std::ofstream output2;
    output.open ("dataP.txt");
    for (int i = 0; i < Nx; ++i) {

        output << p[Nx*Ny*NzHalf + Nx*NyHalf + i] << "\n";
    }
    output.close();

    return 0;
}