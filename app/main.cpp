#include <cmath>
#include "../src/IO/turbulentChannelFlow/writeTurbulentChannelFlowData.h"
#include "../src/interpolation/temporalAdvancement/computeTimeStepOrthogonal.h"
#include "../src/finiteVolumeMethod/fvc/fvc.h"
#include "../src/finiteVolumeMethod/fvm/fvm.h"
#include "../src/finiteVolumeMethod/fvScalarEquation/FvScalarEquation.h"
#include "../src/interpolation/interpolateMDotFromElements2Faces/RhieChowInterpolation.h"
#include "../src/functionObjects/computeBulkVelocity/computeBulkVelocity.h"
#include "../src/functionObjects/computeQCrit/computeQCrit.h"
#include "../src/IO/turbulentChannelFlow/writeTurbulentChannelFlowData2CSV.h"
#include "../src/IO/turbulentChannelFlow/writeTurbulentChannelFlowData2VTK.h"



int main(int argc, char *argv[]) {


    // Mesh parameters
    double delta = 1;                                                   // Reference length
    double Lx = 2*M_PI*delta, Ly = 2*delta, Lz = M_PI*delta;            // Domain size
    int Nx = 32, Ny = 33, Nz = 32;                                      // Number of elements in each direction
    double sx = 0, sy = 3.5, sz = 0;                                    // Hyperbolic tangent mesh stretching


    // Mesh generation
    PolyMesh theMesh;
    theMesh.generatePolyMesh(Lx, Ly, Lz, Nx, Ny, Nz, sx, sy, sz);       // Generate internal mesh
    theMesh.generateBoundaryMesh(Nx, Ny, Nz);                           // Generate boundary mesh


    // Flow properties
    double nu = 1./180;                                                 // Viscosity


    // y-plus of the first node in the wall direction
    double yPlusMin = theMesh.elements[0].centroid.y/nu;


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


    // Turbulent viscosity field initialization (field and BCs)
    ScalarField nut;
    nut.assign(theMesh.nInteriorElements, 0);

    ScalarBoundaryConditions nutBCs;
    nutBCs.addBC("periodic", 0);
    nutBCs.addBC("periodic", 0);
    nutBCs.addBC("fixedValue", 0);
    nutBCs.addBC("fixedValue", 0);
    nutBCs.addBC("periodic", 0);
    nutBCs.addBC("periodic", 0);


    // QCriterion field initialization (field and BCs)
    ScalarField QCrit;
    QCrit.assign(theMesh.nInteriorElements, 0);

    ScalarBoundaryConditions QCritBCs;
    QCritBCs.addBC("zeroGradient", 0);
    QCritBCs.addBC("zeroGradient", 0);
    QCritBCs.addBC("zeroGradient", 0);
    QCritBCs.addBC("zeroGradient", 0);
    QCritBCs.addBC("zeroGradient", 0);
    QCritBCs.addBC("zeroGradient", 0);


    // Linear solver parameters initialization
    LinearSolverConfig pSolver("CGS", 1e-6, 1e6);


    // Transient parameters
    double t = 0;
    double tFinal = 1e24;
    double DeltaT;
    double f = 1;
    double steadyStateCriterion = 1e-4;
    int k = 0, writeInterval = 10000;


    // HPC parameters


    // Pre-definitions
    ScalarField mDot, divUPred, pNew, divUNew;
    VectorField convU, diffU, F, R, RPrev, uPred, gradP, uNew, omega;
    TensorField gradU;
    SparseMatrix laplacianMatrixP;
    FvScalarEquation pEqn;
    double uConvergence, pConvergence, uBulk;
    bool temporalIterate = true;


    // Write the simulation set-up data
    writeTurbulentChannelFlowData(delta, Lx, Ly, Lz, Nx, Ny, Nz, sx, sy, sz, nu, yPlusMin, pSolver.tolerance,
                                  pSolver.solver, pSolver.maxIter, t, tFinal, f, steadyStateCriterion, writeInterval);



    // Time loop FSM algorithm
    while (temporalIterate) {

        // Compute time-step and update time
        DeltaT = computeTimeStepOrthogonal(theMesh, u, nu, f);
        t += DeltaT;
        printf("\nTime = %f s \n", t);


        // Compute convective and diffusive term explicitly
        mDot = RhieChowInterpolation(u, p, theMesh, DeltaT, uBCs, pBCs);
        convU = fvc::convective(mDot, u, theMesh, uBCs);
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
        divUPred = fvc::divergenceNVS(uPred, theMesh, uBCs);
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
        printf("\tThe maximum value of divUNew is: %E \n", divUNew.rms());


        // Compute the difference between the new and old maps and get its maximum-absolute value (to ensure temporal convergence)
        uConvergence = ((uNew - u)/DeltaT).maxAbs();
        pConvergence = ((pNew - p)/DeltaT).maxAbs();
        printf("\tThe steady-state convergence parameter is: %E \n", fmax(uConvergence, pConvergence));


        // Compute the bulk velocity
        uBulk = computeBulkVelocity(uNew, theMesh, uBCs, 0);
        printf("\tThe bulk velocity is: %f \n", uBulk);


        // Compute the vorticity and QCriterion
        gradU = fvc::gradient(u, theMesh, uBCs);
        omega = fvc::curl(gradU, theMesh);
        QCrit = computeQCrit(theMesh, gradU);


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
            writeTurbulentChannelFlowData2CSV(theMesh, pNew, uNew, omega, nut, uBulk, t);
            writeTurbulentChannelFlowData2VTK(theMesh, pNew, uNew, omega, nut, QCrit, pBCs, uBCs, uBCs, nutBCs, QCritBCs, t);
        }

        k++;
    }
    printf("\nSimulation completed!");


    // Write last time-step results to .csv and .VTK file
    printf("\nWriting data corresponding to Time = %f s \n", t);
    writeTurbulentChannelFlowData2CSV(theMesh, pNew, uNew, omega, nut, uBulk, t);
    writeTurbulentChannelFlowData2VTK(theMesh, pNew, uNew, omega, nut, QCrit, pBCs, uBCs, uBCs, nutBCs, QCritBCs, t);


    return 0;
}