//
// Created by ruben on 6/04/24.
//

#include "computeLaminarDHCResults.h"
#include <iostream>
#include <fstream>


void computeLaminarDHCResults(const PolyMesh& theMesh, int Nx, int Ny, int Nz, double Lx, double Ly, double t, const VectorField& u,
                              const ScalarField& T, const ScalarBoundaryConditions& TBCs) {

    // Auxiliary variables definition
    ScalarField Nu;
    double uMaxVerticalMidPlane, y_uMaxVerticalMidPlane;
    double vMaxHorizontalMidPlane, x_vMaxHorizontalMidPlane;
    double NuMean, NuMeanVerticalMidPlane, NuMeanVerticalPlane0;
    double NuMaxVerticalPlane0, y_NuMaxVerticalPlane0;
    double NuMinVerticalPlane0, y_NuMinVerticalPlane0;


    // Computations
    Nu.assign(theMesh.nInteriorElements, 0);
    uMaxVerticalMidPlane = 0, vMaxHorizontalMidPlane = 0;
    NuMean = 0, NuMeanVerticalMidPlane = 0, NuMeanVerticalPlane0 = 0;
    NuMaxVerticalPlane0 = 0, NuMinVerticalPlane0 = 1e24;

    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        // Nusselt number and mean cavity Nusselt number computation
        if ( (theMesh.elements[i].centroid.x != theMesh.elements[0].centroid.x) && (theMesh.elements[i].centroid.x != theMesh.elements[Nx-1].centroid.x) ) {

            Nu[i] = u[i].x*T[i] - (T[i+1] - T[i-1])/(theMesh.elements[i+1].centroid.x - theMesh.elements[i-1].centroid.x);
            NuMean += Nu[i];

        } else if (theMesh.elements[i].centroid.x == theMesh.elements[0].centroid.x) {

            Nu[i] = u[i].x*T[i] - (T[i+1] - TBCs.value[0])/(theMesh.elements[i+1].centroid.x - 0);
            NuMean += Nu[i];

        } else {

            Nu[i] = u[i].x*T[i] - (TBCs.value[1] - T[i-1])/(Lx - theMesh.elements[i-1].centroid.x);
            NuMean += Nu[i];
        }

        // Vertical mid-plane computations
        if (theMesh.elements[i].centroid.x == Lx/2) {

            NuMeanVerticalMidPlane += Nu[i];

            if (u[i].x >= uMaxVerticalMidPlane) {

                uMaxVerticalMidPlane = u[i].x;
                y_uMaxVerticalMidPlane = theMesh.elements[i].centroid.y;
            }
        }

        // Horizontal mid-plane computations
        if (theMesh.elements[i].centroid.y == Ly/2) {

            if (u[i].y >= vMaxHorizontalMidPlane) {

                vMaxHorizontalMidPlane = u[i].y;
                x_vMaxHorizontalMidPlane = theMesh.elements[i].centroid.x;
            }
        }

        // Vertical x=0 plane computations
        if (theMesh.elements[i].centroid.x == theMesh.elements[0].centroid.x) {

            NuMeanVerticalPlane0 += Nu[i];

            if (Nu[i] >= NuMaxVerticalPlane0) {

                NuMaxVerticalPlane0 = Nu[i];
                y_NuMaxVerticalPlane0 = theMesh.elements[i].centroid.y;
            }

            if (Nu[i] <= NuMinVerticalPlane0) {

                NuMinVerticalPlane0 = Nu[i];
                y_NuMinVerticalPlane0 = theMesh.elements[i].centroid.y;
            }
        }
    }
    NuMean /= (Nx*Ny*Nz);
    NuMeanVerticalMidPlane /= (Ny*Nz);
    NuMeanVerticalPlane0 /= (Ny*Nz);


    // Write the computed parameters to .csv file along with time
    std::ofstream outfileData;
    outfileData.open("ResultsDHC.csv", std::ios_base::app);

    if (outfileData.is_open()) {

        // Write the header if file is created
        if (outfileData.tellp() == 0) {
            outfileData << "\"t\"," << "\"uMaxVerticalMidPlane\"," << "\"y_uMaxVerticalMidPlane\","
                        << "\"vMaxHorizontalMidPlane\"," << "\"x_vMaxHorizontalMidPlane\","
                        << "\"NuMean\"," << "\"NuMeanVerticalMidPlane\"," << "\"NuMeanVerticalPlane0\","
                        << "\"NuMaxVerticalPlane0\"," << "\"y_NuMaxVerticalPlane0\","
                        << "\"NuMinVerticalPlane0\"," << "\"y_NuMinVerticalPlane0\"\n";
        }

        // Append data (time, bulk velocity)
        outfileData << std::to_string(t) << "," << uMaxVerticalMidPlane << "," << y_uMaxVerticalMidPlane
                    << "," << vMaxHorizontalMidPlane << "," << x_vMaxHorizontalMidPlane
                    << "," << NuMean << "," << NuMeanVerticalMidPlane << "," << NuMeanVerticalPlane0
                    << "," << NuMaxVerticalPlane0 << "," << y_NuMaxVerticalPlane0
                    << "," << NuMinVerticalPlane0 << "," << y_NuMinVerticalPlane0 << "\n";

    } else {

        printf("Error. Unable to open file. \n");
        std::exit(EXIT_FAILURE);
    }
}