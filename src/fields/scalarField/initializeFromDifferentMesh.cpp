//
// Created by ruben on 10/04/24.
//

#include <string>
#include "initializeFromDifferentMesh.h"
#include "../../IO/in/readCSV.h"


ScalarField interpolateScalarFieldFromDifferentMesh(const PolyMesh& theMesh, double t) {

    // Auxiliary variables initialization
    int iMin1, iMin2;
    double d, dMin1, dMin2, f1, f2;
    Node x1, x2;


    // Initialize the velocity field
    ScalarField T;
    T.assign(theMesh.nInteriorElements, 0);


    // Read the .csv file according to the selected time
    std::vector < std::vector<double> > dataCSV = readCSV("Time_" + std::to_string(t) + ".csv");

    printf("Interpolating VectorField from reference to present mesh... \n");


    // Loop for all the present mesh interior elements
    for (int i = 0; i < theMesh.nInteriorElements; ++i) {

        // Identify the present mesh centroid
        x1 = theMesh.elements[i].centroid;


        // Initialize the minimum distances
        dMin1 = 1e24;
        dMin2 = 1e25;


        // Loop for all the reference mesh interior elements
        for (int j = 0; j < dataCSV.size(); ++j) {

            // Identify the reference mesh centroid
            x2 = {dataCSV[j][0], dataCSV[j][1], dataCSV[j][2]};

            // Compute the distance between both centroids
            d = (x2 - x1).mag();

            // Compare the distance with the minimum ones
            if (d < dMin2) {

                if (d <= dMin1) {

                    dMin2 = dMin1;
                    dMin1 = d;

                    iMin2 = iMin1;
                    iMin1 = j;

                } else {

                    dMin2 = d;
                    iMin2 = j;
                }
            }
        }


        // Assign the closest values from the reference mesh (j) to the present one (i)
        f1 = 1 - dMin1 / (dMin1 + dMin2);
        f2 = 1 - dMin2 / (dMin1 + dMin2);

        T[i] = dataCSV[iMin1][7]*f1 + dataCSV[iMin2][7]*f2;
    }

    printf("Done interpolating! \n");


    return T;
}
