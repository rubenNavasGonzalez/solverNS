//
// Created by ruben on 21/11/23.
//

#include "iterators2ElementVertices.h"


std::vector<int> iterators2ElementVertices(int i, int j, int k, int Nx, int Ny, int Nz) {

    std::vector<int> vertices(8,0);

    vertices[0] = 1*j + (Nx + 1)*i + (Nx + 1)*(Ny + 1)*k;
    vertices[1] = 1*(j + 1) + (Nx + 1)*i + (Nx + 1)*(Ny + 1)*k;
    vertices[2] = 1*(j + 1) + (Nx + 1)*i + (Nx + 1)*(Ny + 1)*(k + 1);
    vertices[3] = 1*j + (Nx + 1)*i + (Nx + 1)*(Ny + 1)*(k + 1);
    vertices[4] = 1*j + (Nx + 1)*(i + 1) + (Nx + 1)*(Ny + 1)*k;
    vertices[5] = 1*(j + 1) + (Nx + 1)*(i + 1) + (Nx + 1)*(Ny + 1)*k;
    vertices[6] = 1*(j + 1) + (Nx + 1)*(i + 1) + (Nx + 1)*(Ny + 1)*(k + 1);
    vertices[7] = 1*j + (Nx + 1)*(i + 1) + (Nx + 1)*(Ny + 1)*(k + 1);

    return vertices;
}