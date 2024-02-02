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


    /*// Mesh parameters
    double delta = 1;                                                       // Reference length
    double Lx = 4*M_PI*delta, Ly = 2*delta, Lz = 4./3*M_PI*delta;           // Domain size
    int m = 8;                                                              // Mesh power factor
    int Nx = pow(2,m), Ny = pow(2,m) + 1, Nz = pow(2,m); // Number of elements in each direction
    double sx = 0, sy = 4.5, sz = 0;                                        // Hyperbolic tangent mesh stretching


    // Mesh generation
    PolyMesh theMesh;
    theMesh.generatePolyMesh(Lx, Ly, Lz, Nx, Ny, Nz, sx, sy, sz);       // Generate internal mesh
    theMesh.generateBoundaryMesh(Nx, Ny, Nz);                           // Generate boundary mesh


    // Velocity field initialization (field and BCs)
    VectorField u;
    u.assign(theMesh.nInteriorElements, {0,0,0});

    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        double x = theMesh.elements[i].centroid.x;
        double y = theMesh.elements[i].centroid.y;
        double z = theMesh.elements[i].centroid.z;
        double pi = M_PI;

        u[i].x = sin(2*pi/Lx*x)*sin(pi/Ly*y)*cos(2*pi/Lz*z);
        u[i].y = -cos(2*pi/Lx*x)*sin(pi/Ly*y)*sin(2*pi/Lz*z);
        u[i].z = sin(2*pi/Lx*x)*sin(pi/Ly*y)*cos(2*pi/Lz*z);
    }

    VectorBoundaryConditions uBCs;
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("fixedValue", {0,0,0});
    uBCs.addBC("fixedValue", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});


    // Compute the vorticity
    TensorField gradU = fvc::gradient(u, theMesh, uBCs);
    VectorField omega = fvc::curl(gradU, theMesh);


    // Verify the gradient field
    ScalarField r11, r12, r13, r21, r22, r23, r31, r32, r33;
    double val11, val12, val13, val21, val22, val23, val31, val32, val33;

    r11.assign(theMesh.nInteriorElements, 0);
    r12.assign(theMesh.nInteriorElements, 0);
    r13.assign(theMesh.nInteriorElements, 0);
    r21.assign(theMesh.nInteriorElements, 0);
    r22.assign(theMesh.nInteriorElements, 0);
    r23.assign(theMesh.nInteriorElements, 0);
    r31.assign(theMesh.nInteriorElements, 0);
    r32.assign(theMesh.nInteriorElements, 0);
    r33.assign(theMesh.nInteriorElements, 0);

    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        double x = theMesh.elements[i].centroid.x;
        double y = theMesh.elements[i].centroid.y;
        double z = theMesh.elements[i].centroid.z;
        double pi = M_PI;

        val11 = (2*pi*cos((2*pi*x)/Lx)*cos((2*pi*z)/Lz)*sin((pi*y)/Ly))/Lx;
        val12 = (2*pi*sin((2*pi*x)/Lx)*sin((pi*y)/Ly)*sin((2*pi*z)/Lz))/Lx;
        val13 = (2*pi*cos((2*pi*x)/Lx)*cos((2*pi*z)/Lz)*sin((pi*y)/Ly))/Lx;
        val21 = (pi*cos((pi*y)/Ly)*cos((2*pi*z)/Lz)*sin((2*pi*x)/Lx))/Ly;
        val22 = -(pi*cos((2*pi*x)/Lx)*cos((pi*y)/Ly)*sin((2*pi*z)/Lz))/Ly;
        val23 = (pi*cos((pi*y)/Ly)*cos((2*pi*z)/Lz)*sin((2*pi*x)/Lx))/Ly;
        val31 = -(2*pi*sin((2*pi*x)/Lx)*sin((pi*y)/Ly)*sin((2*pi*z)/Lz))/Lz;
        val32 = -(2*pi*cos((2*pi*x)/Lx)*cos((2*pi*z)/Lz)*sin((pi*y)/Ly))/Lz;
        val33 = -(2*pi*sin((2*pi*x)/Lx)*sin((pi*y)/Ly)*sin((2*pi*z)/Lz))/Lz;

        r11[i] = fabs(val11 - gradU[i][0][0]);
        r12[i] = fabs(val12 - gradU[i][0][1]);
        r13[i] = fabs(val13 - gradU[i][0][2]);
        r21[i] = fabs(val21 - gradU[i][1][0]);
        r22[i] = fabs(val22 - gradU[i][1][1]);
        r23[i] = fabs(val23 - gradU[i][1][2]);
        r31[i] = fabs(val31 - gradU[i][2][0]);
        r32[i] = fabs(val32 - gradU[i][2][1]);
        r33[i] = fabs(val33 - gradU[i][2][2]);
    }

    printf("\nThe RMS of the (1,1) gradU tensor error is: %f \n", r11.rms());
    printf("The RMS of the (1,2) gradU tensor error is: %f \n", r12.rms());
    printf("The RMS of the (1,3) gradU tensor error is: %f \n", r13.rms());
    printf("The RMS of the (2,1) gradU tensor error is: %f \n", r21.rms());
    printf("The RMS of the (2,2) gradU tensor error is: %f \n", r22.rms());
    printf("The RMS of the (2,3) gradU tensor error is: %f \n", r23.rms());
    printf("The RMS of the (3,1) gradU tensor error is: %f \n", r31.rms());
    printf("The RMS of the (3,2) gradU tensor error is: %f \n", r32.rms());
    printf("The RMS of the (3,3) gradU tensor error is: %f \n", r33.rms());


    // Verify the vorticity field
    ScalarField rX, rY, rZ;
    double valX, valY, valZ;
    rX.assign(theMesh.nInteriorElements, 0);
    rY.assign(theMesh.nInteriorElements, 0);
    rZ.assign(theMesh.nInteriorElements, 0);

    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        double x = theMesh.elements[i].centroid.x;
        double y = theMesh.elements[i].centroid.y;
        double z = theMesh.elements[i].centroid.z;
        double pi = M_PI;

        valX = (pi*cos((2*pi*z)/Lz)*(2*Ly*cos((2*pi*x)/Lx)*sin((pi*y)/Ly) + Lz*cos((pi*y)/Ly)*sin((2*pi*x)/Lx)))/(Ly*Lz);
        valY = -(2*pi*sin((pi*y)/Ly)*(Lz*cos((2*pi*x)/Lx)*cos((2*pi*z)/Lz) + Lx*sin((2*pi*x)/Lx)*sin((2*pi*z)/Lz)))/(Lx*Lz);
        valZ = -(pi*sin((2*pi*x)/Lx)*(Lx*cos((pi*y)/Ly)*cos((2*pi*z)/Lz) - 2*Ly*sin((pi*y)/Ly)*sin((2*pi*z)/Lz)))/(Lx*Ly);

        rX[i] = fabs(valX - omega[i].x);
        rY[i] = fabs(valY - omega[i].y);
        rZ[i] = fabs(valZ - omega[i].z);
    }

    printf("\nThe RMS of the x-vorticity error is: %f \n", rX.rms());
    printf("The RMS of the y-vorticity error is: %f \n", rY.rms());
    printf("The RMS of the z-vorticity error is: %f \n", rZ.rms());*/


    // Mesh parameters
    double delta = 1;                                                   // Reference length
    double Lx = 4*M_PI*delta, Ly = 2*delta, Lz = 4./3*M_PI*delta;       // Domain size
    int Nx = 24, Ny = 25, Nz = 24;                                      // Number of elements in each direction
    double sx = 0, sy = 4, sz = 0;                                      // Hyperbolic tangent mesh stretching


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
    LinearSolverConfig pSolver("CG", 1e-6, 1e6);


    // Transient parameters
    double t = 0;
    double tFinal = 500;
    double DeltaT;
    double f = 1;
    double steadyStateCriterion = 1e-4;
    int k = 0, writeInterval = 2500;


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
        printf("\tThe RMS of divUNew is: %E \n", divUNew.rms());


        // Compute the difference between the new and old maps and get its maximum-absolute value (to ensure temporal convergence)
        uConvergence = ((uNew - u)/DeltaT).maxAbs();
        pConvergence = ((pNew - p)/DeltaT).maxAbs();
        printf("\tThe steady-state convergence parameter is: %E \n", fmax(uConvergence, pConvergence));


        // Compute the bulk velocity
        uBulk = computeBulkVelocity(uNew, theMesh, uBCs, 0);
        printf("\tThe bulk velocity is: %f \n", uBulk);


        // Compute the vorticity and QCriterion
        gradU = fvc::gradient(uNew, theMesh, uBCs);
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
            writeTurbulentChannelFlowData2CSV(theMesh, pNew, uNew, omega, nut, gradU, uBulk, t);
            writeTurbulentChannelFlowData2VTK(theMesh, pNew, uNew, omega, nut, QCrit, pBCs, uBCs, uBCs, nutBCs, QCritBCs, t);
        }

        k++;
    }
    printf("\nSimulation completed!");


    // Write last time-step results to .csv and .VTK file
    printf("\nWriting data corresponding to Time = %f s \n", t);
    writeTurbulentChannelFlowData2CSV(theMesh, pNew, uNew, omega, nut, gradU, uBulk, t);
    writeTurbulentChannelFlowData2VTK(theMesh, pNew, uNew, omega, nut, QCrit, pBCs, uBCs, uBCs, nutBCs, QCritBCs, t);


    return 0;
}