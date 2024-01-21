//
// Created by ruben on 28/12/23.
//

#ifndef SOLVERNS_INITIALIZECHANNELFLOWRETAU180_H
#define SOLVERNS_INITIALIZECHANNELFLOWRETAU180_H


double b = 5.5, kappa = 2.5, A = 0.1*(kappa*log(1/nu) + b);
double uAvg, sinX, sinY, sinZ, cosX, cosY, cosZ;

for (int i = 0; i < theMesh.nInteriorElements; ++i) {

    double y = theMesh.elements[i].centroid.y - 1;

    if ((1 - fabs(y))/nu < 10) {

        uAvg = (1 - fabs(y))/nu;
    } else {

        uAvg = kappa*log((1 - fabs(y))/nu) + b;
    }

    sinX = sin(4*M_PI*theMesh.elements[i].centroid.x/Lx);
    sinY = sin(M_PI*theMesh.elements[i].centroid.y);
    sinZ = sin(2*M_PI*theMesh.elements[i].centroid.z/Lz);
    cosX = cos(4*M_PI*theMesh.elements[i].centroid.x/Lx);
    cosY = 1 + cos(M_PI*theMesh.elements[i].centroid.y);
    cosZ = cos(2*M_PI*theMesh.elements[i].centroid.z/Lz);

    u[i].x = uAvg + A*cosX*sinY*sinZ*Lx/2;
    u[i].y = -A*sinX*cosY*sinZ;
    u[i].z = -A*sinX*sinY*cosZ*Lz/2;
}

#endif //SOLVERNS_INITIALIZECHANNELFLOWRETAU180_H
