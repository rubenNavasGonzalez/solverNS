#include <cmath>
#include "../src/turbulenceModeling/LES/modelLES.h"
#include "../src/IO/out/turbulentChannelFlow/writeTurbulentChannelFlowSetUp.h"
#include "../src/interpolation/temporalAdvancement/computeTimeStepOrthogonal.h"
#include "../src/finiteVolumeMethod/fvc/fvc.h"
#include "../src/finiteVolumeMethod/fvm/fvm.h"
#include "../src/finiteVolumeMethod/fvScalarEquation/FvScalarEquation.h"
#include "../src/interpolation/interpolateMDotFromElements2Faces/RhieChowInterpolation.h"
#include "../src/functionObjects/computeBulkVelocity/computeBulkVelocity.h"
#include "../src/functionObjects/computeQCrit/computeQCrit.h"
#include "../src/IO/out/turbulentChannelFlow/writeTurbulentChannelFlowData2CSV.h"
#include "../src/IO/out/turbulentChannelFlow/writeTurbulentChannelFlowData2VTK.h"


int main() {


    // Mesh parameters
    double delta = 1;                                                           // Reference length
    double Lx = 2*M_PI*delta, Ly = 2*delta, Lz = M_PI*delta;                    // Domain size
    int Nx = 32, Ny = 33, Nz = 32;                                              // Number of elements in each direction
    double sx = 0, sy = 4.5, sz = 0;                                            // Hyperbolic tangent mesh stretching


    // Mesh generation
    PolyMesh theMesh;
    theMesh.generatePolyMesh(Lx, Ly, Lz, Nx, Ny, Nz, sx, sy, sz);               // Generate internal mesh
    theMesh.generateBoundaryMesh(Nx, Ny, Nz);                                   // Generate boundary mesh


    // Flow properties
    double nu = 1./395;                                                         // Viscosity


    // y-plus of the first node in the wall direction
    double yPlusMin = theMesh.elements[0].centroid.y/nu;


    // Linear solver parameters initialization
    LinearSolverConfig pSolver("CG", 1e-6, 1e6);


    // File recording parameters
    int k = 0;                                                                  // Temporal iteration
    double writeIntervalCSV = 1;                                                // Frequency to generate .csv data
    double writeIntervalCSVStatic = writeIntervalCSV;
    double writeIntervalVTK = 10;                                               // Frequency to generate .VTK data
    double writeIntervalVTKStatic = writeIntervalVTK;


    // Transient parameters
    double t = 0;                                                               // Time (dynamic)
    double t0 = t;                                                              // Initial time (static)
    double tFinal = 100;                                                        // Final time
    double DeltaT;                                                              // Time-step
    double f = 1;                                                               // Time-step calculation correction factor
    double steadyStateCriterion = 1e-6;                                         // Steady-state criterion


    // Turbulence modeling
    modelLES turbulenceModel;
    turbulenceModel = None;


    // Velocity field initialization (field and BCs)
    VectorField u;
    u.initialize(theMesh, t, Lx, Ly, Lz, nu, 0);

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


    // Pre-definitions
    ScalarField mDot, divUPred, pNew, divUNew;
    VectorField convU, diffU, F, R, RPrev, uPred, gradP, uNew, omega;
    TensorField gradU;
    SparseMatrix laplacianMatrixP;
    FvScalarEquation pEqn;
    double uConvergence, pConvergence, uBulk;
    bool temporalIterate = true;


    // Write the simulation set-up data
    writeTurbulentChannelFlowSetUp(delta, Lx, Ly, Lz, Nx, Ny, Nz, sx, sy, sz, nu, yPlusMin, pSolver.tolerance,
                                   pSolver.solver, pSolver.maxIter, t, tFinal, f, steadyStateCriterion, writeIntervalCSV,
                                   writeIntervalVTK, turbulenceModel, "caseSetUp");


    // Time loop FSM algorithm
    while (temporalIterate) {


        // Compute the turbulent viscosity
        nut = computeTurbulentViscosity(theMesh, u, uBCs, turbulenceModel);


        // Compute time-step and update time
        DeltaT = computeTimeStepOrthogonal(theMesh, u, nu+nut, f);
        t += DeltaT;
        printf("\nTime = %f s \n", t);


        // Compute convective and diffusive term explicitly
        mDot = RhieChowInterpolation(u, p, theMesh, DeltaT, uBCs, pBCs);
        convU = fvc::convective(mDot, u, theMesh, uBCs);
        diffU = fvc::laplacianOrthogonal(nu+nut, u, theMesh, uBCs);


        // Compute the forcing term   dp/dx = -1   in the x direction
        F = fvc::forcingTerm({-1,0,0}, theMesh);


        // Compute R field and RPrev field (if first time iteration)
        R = diffU - convU - F;

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

        pEqn.constrain(theMesh, pBCs);
        pEqn.perturb();
        pEqn.changeSign();


        // Solve the Poisson Equation with a linear solver to get the new pressure
        printf("\tSolving Poisson Equation...\n");
        pNew = pEqn.solve(pSolver, p);


        // Compute the gradient of the new pressure
        gradP = fvc::gradient(pNew, theMesh, pBCs);


        // Compute the new velocity
        uNew = uPred - (DeltaT/1)*gradP;


        // Compute the divergence of the new velocity explicitly and get its maximum-absolute value (to ensure continuity)
        divUNew = fvc::divergence(uNew, theMesh, uBCs);
        printf("\tThe RMS of divUNew is: %E \n", divUNew.rms());


        // Compute the difference between the new and old maps and get its maximum-absolute value (to ensure temporal convergence)
        uConvergence = ((uNew - u)/DeltaT).maxAbs();
        pConvergence = ((pNew - p)/DeltaT).maxAbs();
        printf("\tThe steady-state convergence parameter is: %E \n", fmax(uConvergence, pConvergence));


        // Compute the bulk velocity
        uBulk = computeBulkVelocity(uNew, theMesh, uBCs, 0, t);
        printf("\tThe bulk velocity is: %f \n", uBulk);


        // Check if new time iteration is needed
        if ( t > tFinal || fmax(uConvergence, pConvergence) < steadyStateCriterion ) {

            temporalIterate = false;
        } else {

            u = uNew;
            p = pNew;
            RPrev = R;
        }


        // Write results to .csv file
        if ( t - t0 >= writeIntervalCSV || k == 0 || !temporalIterate ) {


            // Compute the gradient of velocity tensor and vorticity
            gradU = fvc::gradient(uNew, theMesh, uBCs);
            omega = fvc::curl(gradU, theMesh);

            printf("\nWriting .csv data corresponding to Time = %f s \n", t);
            writeTurbulentChannelFlowData2CSV(theMesh, pNew, uNew, omega, nut, gradU, t);

            if (k != 0) {

                writeIntervalCSV += writeIntervalCSVStatic;
            }
        }


        // Write results to .VTK file
        if ( t - t0 >= writeIntervalVTK || k == 0 || !temporalIterate ) {


            // Compute the gradient of velocity tensor, vorticity and QCriterion
            gradU = fvc::gradient(uNew, theMesh, uBCs);
            omega = fvc::curl(gradU, theMesh);
            QCrit = computeQCrit(theMesh, gradU);

            printf("\nWriting .VTK data corresponding to Time = %f s \n", t);
            writeTurbulentChannelFlowData2VTK(theMesh, pNew, uNew, omega, nut, QCrit, pBCs, uBCs, uBCs, nutBCs, QCritBCs, t);

            if (k != 0) {

                writeIntervalVTK += writeIntervalVTKStatic;
            }
        }

        k++;
    }
    printf("\nSimulation completed! \n");


    return 0;
}