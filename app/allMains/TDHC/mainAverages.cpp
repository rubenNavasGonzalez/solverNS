#include "../src/turbulenceModeling/LES/modelLES.h"
#include "../src/finiteVolumeMethod/fvc/fvc.h"
#include "../src/finiteVolumeMethod/fvm/fvm.h"
#include "../src/finiteVolumeMethod/fvScalarEquation/FvScalarEquation.h"
#include "../src/IO/out/differentiallyHeatedCavity/writeDifferentiallyHeatedCavityData2VTK.h"



int main() {


    // Mesh parameters
    double LRef = 4;                                                            // Reference distance
    double Lx = 1/LRef, Ly = 4/LRef, Lz = 1;                                    // Domain size
    int Nx = 138+1, Ny = 326+1, Nz = 1;                                         // Number of elements in each direction
    double sx = 3.5, sy = 3.5, sz = 0;                                          // Hyperbolic tangent mesh stretching


    // Mesh generation
    PolyMesh theMesh;
    theMesh.generatePolyMesh(Lx, Ly, Lz, Nx, Ny, Nz, sx, sy, sz);               // Generate internal mesh
    theMesh.generateBoundaryMesh(Nx, Ny, Nz);                                   // Generate boundary mesh


    // Transient parameters
    double t = -1;                                                              // Time (dynamic)


    // Velocity field initialization (field and BCs)
    VectorField u;
    u.initialize(theMesh, t, 0);

    VectorBoundaryConditions uBCs;
    uBCs.addBC("fixedValue", {0,0,0});
    uBCs.addBC("fixedValue", {0,0,0});
    uBCs.addBC("fixedValue", {0,0,0});
    uBCs.addBC("fixedValue", {0,0,0});
    uBCs.addBC("empty", {0,0,0});
    uBCs.addBC("empty", {0,0,0});


    // Pressure field initialization (field and BCs)
    ScalarField p;
    p.assign(theMesh.nInteriorElements, 0);

    ScalarBoundaryConditions pBCs;
    pBCs.addBC("zeroGradient", 0);
    pBCs.addBC("zeroGradient", 0);
    pBCs.addBC("zeroGradient", 0);
    pBCs.addBC("zeroGradient", 0);
    pBCs.addBC("empty", 0);
    pBCs.addBC("empty", 0);


    // Temperature field initialization (field and BCs)
    ScalarField T;
    T.initialize(theMesh, t, 0);

    ScalarBoundaryConditions TBCs;
    TBCs.addBC("fixedValue", 1);
    TBCs.addBC("fixedValue", 0);
    TBCs.addBC("zeroGradient", 0);
    TBCs.addBC("zeroGradient", 0);
    TBCs.addBC("empty", 0);
    TBCs.addBC("empty", 0);


    // Turbulent viscosity field initialization (field and BCs)
    ScalarField nut;
    nut.assign(theMesh.nInteriorElements, 0);

    ScalarBoundaryConditions nutBCs;
    nutBCs.addBC("fixedValue", 0);
    nutBCs.addBC("fixedValue", 0);
    nutBCs.addBC("fixedValue", 0);
    nutBCs.addBC("fixedValue", 0);
    nutBCs.addBC("empty", 0);
    nutBCs.addBC("empty", 0);


    // QCriterion field initialization (field and BCs)
    ScalarField QCrit;
    QCrit.assign(theMesh.nInteriorElements, 0);

    ScalarBoundaryConditions QCritBCs;
    QCritBCs.addBC("zeroGradient", 0);
    QCritBCs.addBC("zeroGradient", 0);
    QCritBCs.addBC("zeroGradient", 0);
    QCritBCs.addBC("zeroGradient", 0);
    QCritBCs.addBC("empty", 0);
    QCritBCs.addBC("empty", 0);


    // Pre-definitions
    ScalarField mDot, divUPred, pNew, divUNew, convT, diffT, R_T, RPrev_T, TNew;
    VectorField convU, diffU, buyoancy, R, RPrev, uPred, gradP, uNew, omega;
    omega.assign(theMesh.nInteriorElements, {0,0,0});


    // Write the average field to .VTK format
    writeDifferentiallyHeatedCavityData2VTK(theMesh, p, T, u, omega, nut, QCrit, pBCs, TBCs, uBCs, uBCs, nutBCs, QCritBCs, t);


    return 0;
}