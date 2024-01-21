#include <cmath>
#include "../src/interpolation/temporalAdvancement/computeTimeStepOrthogonal.h"
#include "../src/finiteVolumeMethod/fvc/fvc.h"
#include "../src/finiteVolumeMethod/fvm/fvm.h"
#include "../src/finiteVolumeMethod/fvScalarEquation/FvScalarEquation.h"
#include "../src/interpolation/interpolateMDotFromElements2Faces/RhieChowInterpolation.h"
#include "../src/IO/writeTurbulentChannelFlowData2CSV.h"
#include "../src/IO/writePressureVelocity2VTK.h"
#include "../src/functionObjects/computeBulkVelocity/computeBulkVelocity.h"


int main(int argc, char *argv[]) {

    // Mesh parameters
    double delta = 1;                                                   // Reference length
    double Lx = 4*M_PI*delta, Ly = 2*delta, Lz = 4./3.*M_PI*delta;      // Domain size
    int Nx = 32, Ny = 32, Nz = 32;                                      // Number of elements in each direction
    double sx = 0, sy = 2, sz = 0;                                      // Hyperbolic tangent mesh spacing


    // Mesh generation
    PolyMesh theMesh;
    theMesh.generatePolyMesh(Lx, Ly, Lz, Nx, Ny, Nz, sx, sy, sz);       // Generate internal mesh
    theMesh.generateBoundaryMesh(Nx, Ny, Nz);                           // Generate boundary mesh


    // Flow properties
    double nu = 1./180;                                                  // Viscosity


    // Velocity field initialization (field and BCs)
    VectorField u;
    u.assign(theMesh.nInteriorElements, {0,0,0});
    #include "mainIncludes/initializeChannelFlowReTau180.h"

    VectorBoundaryConditions uBCs;
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("fixedValue", {0,0,0});
    uBCs.addBC("fixedValue", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});


    // Pressure field initialization (field and BCs)
    ScalarField p;
    p.assign(theMesh.nInteriorElements, 0);

    ScalarBoundaryConditions pBCs;
    pBCs.addBC("periodic", 0);
    pBCs.addBC("periodic", 0);
    pBCs.addBC("zeroGradient", 0);
    pBCs.addBC("zeroGradient", 0);
    pBCs.addBC("periodic", 0);
    pBCs.addBC("periodic", 0);


    // Linear solver parameters initialization
    LinearSolverConfig pSolver("BiCGSTAB", 1e-6, 1e6);


    // Transient parameters
    double t = 0;
    double tFinal = 1e24;
    double DeltaT;
    double steadyStateCriterion = 1e-4;
    int k = 0, writeInterval = 5e3;


    // Pre-definitions
    ScalarField mDot, divUPred, pNew, divUNew;
    VectorField convU, diffU, F, R, RPrev, uPred, gradP, uNew, omega;
    TensorField gradU;
    SparseMatrix laplacianMatrixP;
    FvScalarEquation pEqn;
    double uConvergence, pConvergence, uBulk, ReBulk, Cf;
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


        // Compute the forcing term   dp/dx = -1   in the x direction
        F = fvc::forcingTerm({-1,0,0}, theMesh);


        // Compute R field and RPrev field (if first time iteration)
        R = nu*diffU - convU - F;

        if (k == 0) {

            RPrev = R;
        }


        // Compute predictor velocity
        uPred = u + (DeltaT/1)*(1.5*R - 0.5*RPrev);


        // Compute the divergence of UPred explicitly and the laplacian of the pressure field implicitly
        divUPred = fvc::divergence(uPred, theMesh, uBCs);
        laplacianMatrixP = fvm::laplacianOrthogonal(theMesh);


        // Assemble and constrain (apply BCs) the Poisson Equation
        pEqn  =  laplacianMatrixP == (1/DeltaT)*divUPred;
        pEqn.constrain(theMesh, 1/DeltaT, pBCs);


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


        // Compute the gradient of the new velocity field
        gradU = fvc::gradient(uNew, theMesh, uBCs);


        // Compute the vorticity field
        omega = fvc::curl(gradU, theMesh);


        // Compute the difference between the new and old maps and get its maximum-absolute value (to ensure temporal convergence)
        uConvergence = ((uNew - u)/DeltaT).maxAbs();
        pConvergence = ((pNew - p)/DeltaT).maxAbs();
        printf("\tThe steady-state convergence parameter is: %E \n", fmax(uConvergence, pConvergence));


        // Compute the bulk velocity and the bulk Reynolds number
        uBulk = computeBulkVelocity(u, theMesh, uBCs, 0);
        ReBulk = uBulk*(4*Ly*Lz/(2*Ly + 2*Lz))/nu;
        Cf = 1/(0.5*1*pow(uBulk,2));

        printf("\tThe bulk velocity is: %f\n", uBulk);
        printf("\tThe bulk Reynolds number is: %f\n", ReBulk);
        printf("\tThe numerical friction coefficient is: %E\n", Cf);


        // Check if new time iteration is needed
        if ( t > tFinal || fmax(uConvergence, pConvergence) < steadyStateCriterion ) {

            temporalIterate = false;
        } else {

            u = uNew;
            p = pNew;
            RPrev = R;
        }


        // Write results to .csv and .VTK file
        if (k % writeInterval == 0 && k != 0) {

            printf("\nWriting data corresponding to Time = %f s \n", t);
            writeTurbulentChannelFlowData2CSV(pNew, u, omega, uBulk, t);
            writePressureVelocity2VTK(theMesh, pNew, uNew, pBCs, uBCs, t);
        }

        k++;
    }
    printf("\nSimulation completed!");


    // Write last time-step results to .csv and .VTK file
    printf("\nWriting data corresponding to Time = %f s \n", t);
    writeTurbulentChannelFlowData2CSV(pNew, u, omega, uBulk, t);
    writePressureVelocity2VTK(theMesh, pNew, uNew, pBCs, uBCs, t);


    return 0;
}