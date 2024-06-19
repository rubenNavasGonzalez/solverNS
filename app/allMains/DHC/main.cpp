#include <cmath>
#include "../src/turbulenceModeling/LES/modelLES.h"
#include "../src/IO/out/turbulentChannelFlow/writeTurbulentChannelFlowSetUp.h"
#include "../src/interpolation/temporalAdvancement/computeTimeStepOrthogonal.h"
#include "../src/finiteVolumeMethod/fvc/fvc.h"
#include "../src/finiteVolumeMethod/fvm/fvm.h"
#include "../src/finiteVolumeMethod/fvScalarEquation/FvScalarEquation.h"
#include "../src/interpolation/interpolateMDotFromElements2Faces/RhieChowInterpolation.h"
#include "../src/functionObjects/computeQCrit/computeQCrit.h"
#include "../src/IO/out/differentiallyHeatedCavity/writeDifferentiallyHeatedCavityData2VTK.h"
#include "../src/functionObjects/computeLaminarDHCResults/computeLaminarDHCResults.h"


int main() {


    // Mesh parameters
    double delta = 1;                                                           // Reference length
    double Lx = 1*delta, Ly = 1*delta, Lz = 1*delta;                            // Domain size
    int Nx = 128+1, Ny = 128+1, Nz = 1;                                         // Number of elements in each direction
    double sx = 0, sy = 0, sz = 0;                                              // Hyperbolic tangent mesh stretching


    // Mesh generation
    PolyMesh theMesh;
    theMesh.generatePolyMesh(Lx, Ly, Lz, Nx, Ny, Nz, sx, sy, sz);               // Generate internal mesh
    theMesh.generateBoundaryMesh(Nx, Ny, Nz);                                   // Generate boundary mesh


    // Flow properties
    double Pr = 0.71;                                                           // Prandtl number
    double Ra = 1e6;                                                            // Rayleigh number


    // Linear solver parameters initialization
    LinearSolverConfig pSolver("CG", 1e-6, 1e6);


    // File recording parameters
    int k = 0;                                                                  // Temporal iteration
    double writeIntervalCSV = 0.1;                                              // Frequency to generate .csv data
    double writeIntervalCSVStatic = writeIntervalCSV;
    double writeIntervalVTK = 1e24;                                             // Frequency to generate .VTK data
    double writeIntervalVTKStatic = writeIntervalVTK;


    // Transient parameters
    double t = 0;                                                               // Time (dynamic)
    double t0 = t;                                                              // Initial time (static)
    double tFinal = 1e24;                                                       // Final time
    double DeltaT;                                                              // Time-step
    double f = 0.5;                                                             // Time-step calculation correction factor
    double steadyStateCriterion = 1e-24;                                        // Steady-state criterion


    // Turbulence modeling
    modelLES turbulenceModel;
    turbulenceModel = None;


    // Velocity field initialization (field and BCs)
    VectorField u;
    u.assign(theMesh.nInteriorElements, {0,0,0});

    VectorBoundaryConditions uBCs;
    uBCs.addBC("fixedValue", {0,0,0});
    uBCs.addBC("fixedValue", {0,0,0});
    uBCs.addBC("fixedValue", {0,0,0});
    uBCs.addBC("fixedValue", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});


    // Pressure field initialization (field and BCs)
    ScalarField p;
    p.assign(theMesh.nInteriorElements, 0);

    ScalarBoundaryConditions pBCs;
    pBCs.addBC("zeroGradient", 0);
    pBCs.addBC("zeroGradient", 0);
    pBCs.addBC("zeroGradient", 0);
    pBCs.addBC("zeroGradient", 0);
    pBCs.addBC("periodic", 0);
    pBCs.addBC("periodic", 0);


    // Temperature field initialization (field and BCs)
    ScalarField T;
    T.assign(theMesh.nInteriorElements, 0);

    ScalarBoundaryConditions TBCs;
    TBCs.addBC("fixedValue", 1);
    TBCs.addBC("fixedValue", 0);
    TBCs.addBC("zeroGradient", 0);
    TBCs.addBC("zeroGradient", 0);
    TBCs.addBC("periodic", 0);
    TBCs.addBC("periodic", 0);


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
    ScalarField mDot, divUPred, pNew, divUNew, convT, diffT, R_T, RPrev_T, TNew;
    VectorField convU, diffU, buyoancy, R, RPrev, uPred, gradP, uNew, omega;
    TensorField gradU;
    SparseMatrix laplacianMatrixP;
    FvScalarEquation pEqn;
    double uConvergence, pConvergence, TConvergence;
    bool temporalIterate = true;


    // Write the simulation set-up data
    /*writeTurbulentChannelFlowSetUp(delta, Lx, Ly, Lz, Nx, Ny, Nz, sx, sy, sz, nu, yPlusMin, pSolver.tolerance,
                                   pSolver.solver, pSolver.maxIter, t, tFinal, f, steadyStateCriterion, writeIntervalCSV,
                                   writeIntervalVTK, turbulenceModel, "caseSetUp");*/


    // Time loop FSM algorithm
    while (temporalIterate) {


        // Compute the turbulent viscosity
        nut = computeTurbulentViscosity(theMesh, u, uBCs, turbulenceModel);


        // Compute time-step and update time
        DeltaT = computeTimeStepOrthogonal(theMesh, u, Pr, f);
        t += DeltaT;
        printf("\nTime = %f s \n", t);


        // Compute convective and diffusive term explicitly
        mDot = RhieChowInterpolation(u, p, theMesh, DeltaT, uBCs, pBCs);
        convU = fvc::convective(mDot, u, theMesh, uBCs);
        diffU = fvc::laplacianOrthogonal(u, theMesh, uBCs);


        // Compute the buyoancy term
        buyoancy = fvc::buyoancy(theMesh, Pr*Ra, T, {0,1,0});


        // Compute R field and RPrev field (if first time iteration)
        R = Pr*diffU - convU + buyoancy;

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


        // Compute convective and diffusive term of the temperature equation
        convT = fvc::convective(mDot, T, theMesh, TBCs);
        diffT = fvc::laplacianOrthogonal(T, theMesh, TBCs);


        // Compute R_T field and RPrev_T field (if first time iteration)
        R_T = diffT - convT;

        if (k == 0) {

            RPrev_T = R_T;
        }


        // Compute the new temperature
        TNew = T + (DeltaT/1)*(1.5*R_T - 0.5*RPrev_T);


        // Compute the difference between the new and old maps and get its maximum-absolute value (to ensure temporal convergence)
        uConvergence = ((uNew - u)/DeltaT).maxAbs();
        pConvergence = ((pNew - p)/DeltaT).maxAbs();
        TConvergence = ((TNew - T)/DeltaT).maxAbs();
        printf("\tThe steady-state convergence parameter is: %E \n", fmax(uConvergence, fmax(pConvergence, TConvergence)));


        // Check if new time iteration is needed
        if ( t > tFinal || fmax(uConvergence, fmax(pConvergence, TConvergence)) < steadyStateCriterion ) {

            temporalIterate = false;
        } else {

            u = uNew;
            p = pNew;
            T = TNew;
            RPrev = R;
            RPrev_T = R_T;
        }


        // Write results to .csv file
        if ( t - t0 >= writeIntervalCSV || k == 0 || !temporalIterate ) {

            printf("\nWriting .csv data corresponding to Time = %f s \n", t);
            computeLaminarDHCResults(theMesh, Nx, Ny, Nz, Lx, Ly, t, u, T, TBCs);

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
            writeTurbulentChannelFlowData2VTK(theMesh, p, T, u, omega, nut, QCrit, pBCs, TBCs, uBCs, uBCs, nutBCs, QCritBCs, t);

            if (k != 0) {

                writeIntervalVTK += writeIntervalVTKStatic;
            }
        }

        k++;
    }
    printf("\nSimulation completed! \n");


    return 0;
}