#include <cmath>
#include "../src/interpolation/temporalAdvancement/computeTimeStepOrthogonal.h"
#include "../src/finiteVolumeMethod/fvc/fvc.h"
#include "../src/finiteVolumeMethod/fvm/fvm.h"
#include "../src/finiteVolumeMethod/fvScalarEquation/FvScalarEquation.h"
#include "../src/interpolation/interpolateMDotFromElements2Faces/RhieChowInterpolation.h"


int main(int argc, char *argv[]) {

    // Mesh parameters
    double Lx = 10, Ly = 1, Lz = 1;                                     // Domain size
    int Nx = 128, Ny = 32, Nz = 8;                                      // Number of elements in each direction
    double sx = 0, sy = 2.5, sz = 0;                                    // Hyperbolic tangent mesh spacing


    // Mesh generation
    PolyMesh theMesh;
    theMesh.generatePolyMesh(Lx, Ly, Lz, Nx, Ny, Nz, sx, sy, sz);       // Generate internal mesh
    theMesh.generateBoundaryMesh(Nx, Ny, Nz);                           // Generate boundary mesh


    // Flow properties
    double nu = 0.05;                                                   // Viscosity


    // Velocity field initialization (field and BCs)
    VectorField u;
    u.assign(theMesh.nInteriorElements, {0,0,0});

    VectorBoundaryConditions uBCs;
    uBCs.addBC("fixedValue", {1,0,0});
    uBCs.addBC("zeroGradient", {0,0,0});
    uBCs.addBC("fixedValue", {0,0,0});
    uBCs.addBC("fixedValue", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});


    // Pressure field initialization (field and BCs)
    ScalarField p;
    p.assign(theMesh.nInteriorElements, 0);

    ScalarBoundaryConditions pBCs;
    pBCs.addBC("zeroGradient", 0);
    pBCs.addBC("fixedValue", 0);
    pBCs.addBC("zeroGradient", 0);
    pBCs.addBC("zeroGradient", 0);
    pBCs.addBC("periodic", 0);
    pBCs.addBC("periodic", 0);


    // Linear solver parameters initialization
    LinearSolverConfig pSolver("BiCGSTAB", 1e-6, 1e6);


    // Transient parameters
    double t = 0;
    double tFinal = 10;
    double DeltaT;
    double steadyStateCriterion = 1e-3;


    // Pre-definitions
    ScalarField mDot, divUPred, pNew, divUNew;
    VectorField convU, diffU, R, RPrev, uPred, gradP, uNew;
    SparseMatrix laplacianMatrixP;
    FvScalarEquation pEqn;
    double uConvergence, pConvergence;
    bool temporalIterate = true;



    // Time loop FSM algorithm
    while (temporalIterate) {

        // Compute time-step and update time
        DeltaT = computeTimeStepOrthogonal(theMesh, u, nu);
        t += DeltaT;
        printf("\nTime = %f s \n", t);


        // Compute convective and diffusive term explicitly
        mDot = RhieChowInterpolation(u, p, theMesh, DeltaT, uBCs, pBCs);
        convU = fvc::convectiveOrthogonal(mDot, u, theMesh, uBCs, "CDS");
        diffU = fvc::laplacianOrthogonal(u, theMesh, uBCs);


        // Compute R field and RPrev field (if first time iteration)
        R = nu*diffU - convU;

        if (t == DeltaT) {

            RPrev = R;
        }


        // Compute predictor velocity
        uPred = u + (DeltaT/1)*(1.5*R - 0.5*RPrev);


        // Compute the divergence of UPred explicitly and the laplacian of the pressure field implicitly
        divUPred = fvc::divergence(uPred, theMesh, uBCs);
        laplacianMatrixP = fvm::laplacianOrthogonal(theMesh);


        // Assemble and constrain (apply BCs) the Poisson Equation
        pEqn  =  laplacianMatrixP == (1/DeltaT)*divUPred;
        pEqn.constrain(theMesh, pBCs);


        // Solve the Poisson Equation with a linear solver to get the new pressure
        printf("\tSolving Poisson Equation...\n");
        pNew = pEqn.solve(pSolver, p);


        // Compute the gradient of the new pressure
        gradP = fvc::gradient(pNew, theMesh, pBCs);


        // Compute the new velocity
        uNew = uPred - (DeltaT/1)*gradP;


        // Compute the divergence of the new velocity explicitly and get its maximum-absolute value (to ensure continuity)
        divUNew = fvc::divergence(uNew, theMesh, uBCs);
        printf("\tThe maximum value of divUNew is: %E \n", divUNew.max());


        // Compute the difference between the new and old maps and get its maximum-absolute value (to ensure temporal convergence)
        uConvergence = ((uNew - u)/DeltaT).maxAbs();
        pConvergence = ((pNew - p)/DeltaT).maxAbs();
        printf("\tThe steady-state convergence parameter is: %E \n", fmax(uConvergence, pConvergence));


        // Check if new time iteration is needed
        if (t > tFinal || fmax(uConvergence, pConvergence) < steadyStateCriterion) {

            temporalIterate = false;
        } else {

            u = uNew;
            p = pNew;
            RPrev = R;
        }
    }


    // Write velocity and pressure fields to a .vtk file
    u.writeVectorField2VTK("U", theMesh, uBCs);
    p.writeScalarField2VTK("p", theMesh, pBCs);

    printf("Simulation completed!");


    return 0;
}