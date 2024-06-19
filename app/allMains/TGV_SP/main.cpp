#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "../src/interpolation/temporalAdvancement/computeTimeStepOrthogonal.h"
#include "../src/finiteVolumeMethod/fvc/fvc.h"
#include "../src/finiteVolumeMethod/fvm/fvm.h"
#include "../src/finiteVolumeMethod/fvScalarEquation/FvScalarEquation.h"
#include "../src/interpolation/interpolateMDotFromElements2Faces/RhieChowInterpolation.h"



int main() {


    // Mesh parameters
    double L = 1;                                                               // Reference length
    double Lx = 2*M_PI*L, Ly = 2*M_PI*L, Lz = 2*M_PI*L;                         // Domain size
    int Nx = 64, Ny = 64, Nz = 64;                                              // Number of elements in each direction
    double sx = 0, sy = 0, sz = 0;                                              // Hyperbolic tangent mesh stretching
    

    // Mesh generation
    PolyMesh theMesh;
    theMesh.generatePolyMesh(Lx, Ly, Lz, Nx, Ny, Nz, sx, sy, sz);               // Generate internal mesh
    theMesh.generateBoundaryMesh(Nx, Ny, Nz);                                   // Generate boundary mesh


    // Flow properties
    double nu = 1./1600;                                                        // Viscosity


    // Linear solver parameters initialization
    LinearSolverConfig pSolver("CG", 1e-9, 1e6);


    // Transient parameters
    double t = 0;                                                               // Time (dynamic)
    double tFinal = 20;                                                         // Final time
    double DeltaT;                                                              // Time-step
    double f = 1;                                                               // Time-step factor
    double steadyStateCriterion = 1e-6;                                         // Steady-state criterion
    int k = 0;                                                                  // Temporal iteration


    // Velocity field initialization (field and BCs)
    VectorField u, v;
    u.initialize(theMesh, t, Lx, Ly, Lz, 0);

    VectorBoundaryConditions uBCs;
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});
    uBCs.addBC("periodic", {0,0,0});

    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        v.push_back({double(i), 2*double(i), 3*double(i)});
    }


    // Pressure field initialization (field and BCs)
    ScalarField p;
    p.assign(theMesh.nInteriorElements, 0);

    ScalarBoundaryConditions pBCs;
    pBCs.addBC("periodic", 0);
    pBCs.addBC("periodic", 0);
    pBCs.addBC("periodic", 0);
    pBCs.addBC("periodic", 0);
    pBCs.addBC("periodic", 0);
    pBCs.addBC("periodic", 0);


    // Turbulent viscosity field initialization (field and BCs)
    ScalarField nut;
    nut.assign(theMesh.nInteriorElements, 0);


    // Pre-definitions
    ScalarField mDot, divUPred, pNew, divUNew, divU, divUCorrected;
    VectorField convU, diffU, diffV, diffU2, R, RPrev, uPred, gradP, uNew;
    TensorField gradU;
    SparseMatrix laplacianMatrixP;
    FvScalarEquation pEqn;
    double convSkewSymmetry, diffSymmetric, diffPositiveDefined, divUMag, divUCorrectedMag;
    double Ek_n, Ek_n1, RHS, derEk_conv, derEk_diff, derEk_grad;
    double uConvergence, pConvergence;
    bool temporalIterate = true;



    // Time loop FSM algorithm
    while (temporalIterate) {


        // Update time
        //DeltaT = computeTimeStepOrthogonal(theMesh, u, nu + nut, f);
        DeltaT = 0.025;
        t += DeltaT;
        printf("\nTime = %f s \n", t);


        // Compute convective and diffusive term explicitly
        mDot = RhieChowInterpolation(u, p, theMesh, DeltaT, uBCs, pBCs);
        convU = fvc::convective(mDot, u, theMesh, uBCs);
        diffU = fvc::laplacianOrthogonal(nu+nut, u, theMesh, uBCs);
        diffV = fvc::laplacianOrthogonal(nu+nut, v, theMesh, uBCs);


        // Differential operators symmetries
        convSkewSymmetry = 0;
        diffSymmetric = 0;
        diffPositiveDefined = 0;

        for (int i = 0; i < theMesh.nInteriorElements; ++i) {

            convSkewSymmetry += u[i].x * theMesh.elements[i].Vf * convU[i].x;
            convSkewSymmetry += u[i].y * theMesh.elements[i].Vf * convU[i].y;
            convSkewSymmetry += u[i].z * theMesh.elements[i].Vf * convU[i].z;

            diffSymmetric += v[i].x * theMesh.elements[i].Vf * diffU[i].x - u[i].x * theMesh.elements[i].Vf * diffV[i].x;
            diffSymmetric += v[i].y * theMesh.elements[i].Vf * diffU[i].y - u[i].y * theMesh.elements[i].Vf * diffV[i].y;
            diffSymmetric += v[i].z * theMesh.elements[i].Vf * diffU[i].z - u[i].z * theMesh.elements[i].Vf * diffV[i].z;

            diffPositiveDefined += u[i].x * theMesh.elements[i].Vf * diffU[i].x;
            diffPositiveDefined += u[i].y * theMesh.elements[i].Vf * diffU[i].y;
            diffPositiveDefined += u[i].z * theMesh.elements[i].Vf * diffU[i].z;
        }
        convSkewSymmetry = fabs(convSkewSymmetry);
        diffSymmetric = fabs(diffSymmetric);

        divU = fvc::divergence(u, theMesh, uBCs);
        divUMag = divU.rms();

        divUCorrected.assign(theMesh.nElements, 0);

        for (int i = 0; i < theMesh.nInteriorFaces; ++i) {

            divUCorrected[ theMesh.faces[i].iOwner ] += mDot[i];
            divUCorrected[ theMesh.faces[i].iNeighbour ] -= mDot[i];
        }
        for (int j = 0; j < theMesh.nBoundaries; ++j) {


            for (int i = theMesh.boundaries[j].startFace; i < theMesh.boundaries[j].startFace + theMesh.boundaries[j].nBoundaryFaces; ++i) {

                divUCorrected[ theMesh.faces[i].iOwner ] += mDot[i];
            }
        }

        divUCorrectedMag = divUCorrected.rms();


        // Write differential operators symmetries to .csv file
        std::ofstream outfileSymmetry;
        outfileSymmetry.open("Symmetry.csv", std::ios_base::app);

        if (outfileSymmetry.tellp() == 0) {

            outfileSymmetry << "\"t\"," << "\"convSkewSymmetric\"," << "\"diffSymmetric\"," << "\"diffPositiveDefined\","
                            << "\"divUMag\"," << "\"divUCorrectedMag\"\n";
        }
        outfileSymmetry << t - DeltaT << "," << convSkewSymmetry << "," << diffSymmetric << "," << diffPositiveDefined
                        << "," << divUMag << "," << divUCorrectedMag << "\n";
        outfileSymmetry.close();


        // Compute R field and RPrev field (if first time iteration)
        R = diffU - convU;

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


        // Momentum equation energetic analysis
        Ek_n = 0;
        Ek_n1 = 0;
        RHS = 0;
        derEk_conv = 0;
        derEk_diff = 0;
        derEk_grad = 0;

        gradP = fvc::gradient(p, theMesh, pBCs);

        for (int i = 0; i < theMesh.nInteriorElements; ++i) {

            Ek_n += 0.5 * theMesh.elements[i].Vf * (pow(u[i].x ,2) + pow(u[i].y ,2) + pow(u[i].z ,2));
            Ek_n1 += 0.5 * theMesh.elements[i].Vf * (pow(uNew[i].x ,2) + pow(uNew[i].y ,2) + pow(uNew[i].z ,2));

            RHS += DeltaT * 0.5 * (u[i].x + uNew[i].x) * theMesh.elements[i].Vf * (R[i].x - gradP[i].x);
            RHS += DeltaT * 0.5 * (u[i].y + uNew[i].y) * theMesh.elements[i].Vf * (R[i].y - gradP[i].y);
            RHS += DeltaT * 0.5 * (u[i].z + uNew[i].z) * theMesh.elements[i].Vf * (R[i].z - gradP[i].z);

            derEk_conv -= 0.5 * (u[i].x + uNew[i].x) * theMesh.elements[i].Vf * convU[i].x;
            derEk_conv -= 0.5 * (u[i].y + uNew[i].y) * theMesh.elements[i].Vf * convU[i].y;
            derEk_conv -= 0.5 * (u[i].z + uNew[i].z) * theMesh.elements[i].Vf * convU[i].z;

            derEk_diff += 0.5 * (u[i].x + uNew[i].x) * theMesh.elements[i].Vf * diffU[i].x;
            derEk_diff += 0.5 * (u[i].y + uNew[i].y) * theMesh.elements[i].Vf * diffU[i].y;
            derEk_diff += 0.5 * (u[i].z + uNew[i].z) * theMesh.elements[i].Vf * diffU[i].z;

            derEk_grad -= 0.5 * (u[i].x + uNew[i].x) * theMesh.elements[i].Vf * gradP[i].x;
            derEk_grad -= 0.5 * (u[i].y + uNew[i].y) * theMesh.elements[i].Vf * gradP[i].y;
            derEk_grad -= 0.5 * (u[i].z + uNew[i].z) * theMesh.elements[i].Vf * gradP[i].z;
        }

        // Write the momentum equation energetic analysis to a file
        std::ofstream outfileEnergy;
        outfileEnergy.open("Energy.csv", std::ios_base::app);

        if (outfileEnergy.tellp() == 0) {

            outfileEnergy << "\"t_n\"," << "\"t_n1\"," << "\"DeltaT\"," << "\"Ek_n\"," << "\"Ek_n1\"," << "\"RHS\","
            << "\"derEk_conv\"," << "\"derEk_diff\"," << "\"derEk_grad\"\n";
        }
        outfileEnergy << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1)
        << t - DeltaT << "," << t << "," << DeltaT << "," << Ek_n << "," << Ek_n1 << "," << RHS << ","
        << derEk_conv << "," << derEk_diff << "," << derEk_grad << "\n";
        outfileEnergy.close();

        printf("\tThe energy residual is: %E \n", fabs(Ek_n1 - Ek_n - RHS));


        // Compute the divergence of the new velocity explicitly and get its maximum-absolute value (to ensure continuity)
        divUNew = fvc::divergence(uNew, theMesh, uBCs);
        printf("\tThe RMS of divUNew is: %E \n", divUNew.rms());


        // Compute the difference between the new and old maps and get its maximum-absolute value (to ensure temporal convergence)
        uConvergence = ((uNew - u)/DeltaT).maxAbs();
        pConvergence = ((pNew - p)/DeltaT).maxAbs();
        printf("\tThe steady-state convergence parameter is: %E \n", fmax(uConvergence, pConvergence));


        // Check if new time iteration is needed
        if ( t > tFinal || fmax(uConvergence, pConvergence) < steadyStateCriterion ) {

            temporalIterate = false;
        } else {

            u = uNew;
            p = pNew;
            RPrev = R;
        }
        k++;
    }
    printf("\nSimulation completed! \n");


    return 0;
}