//
// Created by ruben on 23/11/23.
//

#include <cstdio>
#include "PolyMesh.h"
#include "iterators2ElementVertices.h"
#include "computeRectangleArea.h"


void PolyMesh::generateBoundaryMesh(int Nx, int Ny, int Nz) {

    printf("Generating boundary mesh...\n");


    // Auxiliary variables declaration
    int startFaceRight, startFaceLeft, startFaceBottom, startFaceTop, startFaceBack, startFaceFront;
    std::vector<int> vertices;
    Element auxElement;
    Face auxFace;
    Boundary auxBoundary;


    // Computation of the starting face of each boundary
    startFaceLeft = nInteriorFaces;
    startFaceRight = startFaceLeft + Ny*Nz;
    startFaceBottom = startFaceRight + Ny*Nz;
    startFaceTop = startFaceBottom + Nx*Nz;
    startFaceBack = startFaceTop + Nx*Nz;
    startFaceFront = startFaceBack + Nx*Ny;


    // Generation of the x-perpendicular left boundary
    auxBoundary.startFace = startFaceLeft;
    auxBoundary.nBoundaryFaces = Ny*Nz;
    boundaries.push_back(auxBoundary);
    nBoundaries++;

    for (int k = 0; k < Nz; ++k) {
        for (int i = 0; i < Ny; ++i) {

            vertices = iterators2ElementVertices(i, 0, k, Nx, Ny, Nz);

            // Computation of the vertices of each face
            auxFace.iNodes = {vertices[0], vertices[4], vertices[7], vertices[3]};

            // Computation of the owner and the neighbour element (face to element connectivity) and periodic face
            // It is assumed that the owner is the element with the lowest index
            auxFace.iOwner = Nx*i + Nx*Ny*k;
            auxFace.iNeighbour = nElements;
            faces[(Nx - 1)*i + (Nx - 1)*Ny*k].iOwnerFar = nElements;
            faces[(Nx - 1)*i + (Nx - 1)*Ny*k].iNeighbourFar = faces[(Nx - 1)*i + (Nx - 1)*Ny*k].iNeighbour + 1;
            auxFace.iPeriodicFace = nFaces + auxBoundary.nBoundaryFaces;

            // Computation of the surface and the surface vector. The surface vector points to the neighbour element
            auxFace.SfMag = computeRectangleArea(nodes[vertices[0]].y, nodes[vertices[4]].y,
                                                 nodes[vertices[0]].z, nodes[vertices[3]].z);
            auxFace.Sf = {-auxFace.SfMag, 0, 0};

            // Add the auxiliary face to the mesh faces array
            faces.push_back(auxFace);
            nBoundaryFaces++;
            nFaces++;

            // Computation of the vertices of each element
            auxElement.iNodes = {vertices[0], vertices[4], vertices[7], vertices[3]};

            // Add the auxiliary element to the mesh vector of elements
            elements.push_back(auxElement);
            nBoundaryElements++;
            nElements++;
        }
    }


    // Generation of the x-perpendicular right boundary
    auxBoundary.startFace = startFaceRight;
    auxBoundary.nBoundaryFaces = Ny*Nz;
    boundaries.push_back(auxBoundary);
    nBoundaries++;

    for (int k = 0; k < Nz; ++k) {
        for (int i = 0; i < Ny; ++i) {

            vertices = iterators2ElementVertices(i, Nx - 1, k, Nx, Ny, Nz);

            // Computation of the vertices of each face
            auxFace.iNodes = {vertices[1], vertices[5], vertices[6], vertices[2]};

            // Computation of the owner and the neighbour element (face to element connectivity) and periodic face
            // It is assumed that the owner is the element with the lowest index
            auxFace.iOwner = (Nx - 1) + Nx*i + Nx*Ny*k;
            auxFace.iNeighbour = nElements;
            faces[Nx - 2 + (Nx - 1)*i + (Nx - 1)*Ny*k].iOwnerFar = faces[Nx - 2 + (Nx - 1)*i + (Nx - 1)*Ny*k].iOwner - 1;
            faces[Nx - 2 + (Nx - 1)*i + (Nx - 1)*Ny*k].iNeighbourFar = nElements;
            auxFace.iPeriodicFace = nFaces - auxBoundary.nBoundaryFaces;

            // Computation of the surface and the surface vector. The surface vector points to the neighbour element
            auxFace.SfMag = computeRectangleArea(nodes[vertices[1]].y, nodes[vertices[5]].y,
                                                 nodes[vertices[1]].z, nodes[vertices[2]].z);
            auxFace.Sf = {auxFace.SfMag, 0, 0};

            // Add the auxiliary face to the mesh faces array
            faces.push_back(auxFace);
            nBoundaryFaces++;
            nFaces++;

            // Computation of the vertices of each element
            auxElement.iNodes = {vertices[1], vertices[5], vertices[6], vertices[2]};

            // Add the auxiliary element to the mesh vector of elements
            elements.push_back(auxElement);
            nBoundaryElements++;
            nElements++;
        }
    }


    // Generation of the y-perpendicular bottom boundary
    auxBoundary.startFace = startFaceBottom;
    auxBoundary.nBoundaryFaces = Nx*Nz;
    boundaries.push_back(auxBoundary);
    nBoundaries++;

    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Nx; ++j) {

            vertices = iterators2ElementVertices(0, j, k, Nx, Ny, Nz);

            // Computation of the vertices of each face
            auxFace.iNodes = {vertices[0], vertices[1], vertices[2], vertices[3]};

            // Computation of the owner and the neighbour element (face to element connectivity) and periodic face
            // It is assumed that the owner is the element with the lowest index
            auxFace.iOwner = j + Nx*Ny*k;
            auxFace.iNeighbour = nElements;
            faces[(Nx - 1)*Ny*Nz + j + Nx*(Ny - 1)*k].iOwnerFar = nElements;
            faces[(Nx - 1)*Ny*Nz + j + Nx*(Ny - 1)*k].iNeighbourFar = faces[(Nx - 1)*Ny*Nz + j + Nx*(Ny - 1)*k].iNeighbour + Nx;
            auxFace.iPeriodicFace = nFaces + auxBoundary.nBoundaryFaces;

            // Computation of the surface and the surface vector. The surface vector points to the neighbour element
            auxFace.SfMag = computeRectangleArea(nodes[vertices[0]].x, nodes[vertices[1]].x,
                                                 nodes[vertices[0]].z, nodes[vertices[3]].z);
            auxFace.Sf = {0, -auxFace.SfMag, 0};

            // Add the auxiliary face to the mesh faces array
            faces.push_back(auxFace);
            nBoundaryFaces++;
            nFaces++;

            // Computation of the vertices of each element
            auxElement.iNodes = {vertices[0], vertices[1], vertices[2], vertices[3]};

            // Add the auxiliary element to the mesh vector of elements
            elements.push_back(auxElement);
            nBoundaryElements++;
            nElements++;
        }
    }


    // Generation of the y-perpendicular top boundary
    auxBoundary.startFace = startFaceTop;
    auxBoundary.nBoundaryFaces = Nx*Nz;
    boundaries.push_back(auxBoundary);
    nBoundaries++;

    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Nx; ++j) {

            vertices = iterators2ElementVertices(Ny - 1, j, k, Nx, Ny, Nz);

            // Computation of the vertices of each face
            auxFace.iNodes = {vertices[4], vertices[5], vertices[6], vertices[7]};

            // Computation of the owner and the neighbour element (face to element connectivity) and periodic face
            // It is assumed that the owner is the element with the lowest index
            auxFace.iOwner = Nx*(Ny - 1) + j + Nx*Ny*k;
            auxFace.iNeighbour = nElements;
            faces[(Nx - 1)*Ny*Nz + Nx*(Ny - 2) + j + Nx*(Ny - 1)*k].iOwnerFar = faces[(Nx - 1)*Ny*Nz + Nx*(Ny - 2) + j + Nx*(Ny - 1)*k].iOwner - Nx;
            faces[(Nx - 1)*Ny*Nz + Nx*(Ny - 2) + j + Nx*(Ny - 1)*k].iNeighbourFar = nElements;
            auxFace.iPeriodicFace = nFaces - auxBoundary.nBoundaryFaces;

            // Computation of the surface and the surface vector. The surface vector points to the neighbour element
            auxFace.SfMag = computeRectangleArea(nodes[vertices[4]].x, nodes[vertices[5]].x,
                                                 nodes[vertices[4]].z, nodes[vertices[7]].z);
            auxFace.Sf = {0, auxFace.SfMag, 0};

            // Add the auxiliary face to the mesh faces array
            faces.push_back(auxFace);
            nBoundaryFaces++;
            nFaces++;

            // Computation of the vertices of each element
            auxElement.iNodes = {vertices[4], vertices[5], vertices[6], vertices[7]};

            // Add the auxiliary element to the mesh vector of elements
            elements.push_back(auxElement);
            nBoundaryElements++;
            nElements++;
        }
    }


    printf("Boundary mesh generated successfully!!\n");
}