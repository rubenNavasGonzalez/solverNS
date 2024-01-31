//
// Created by ruben on 19/11/23.
//

#include <cstdio>
#include <cmath>
#include "PolyMesh.h"
#include "iterators2ElementVertices.h"
#include "computeHexahedronVolume.h"
#include "computeHexahedronCentroid.h"
#include "computeRectangleCentroid.h"
#include "computeRectangleArea.h"


void PolyMesh::generatePolyMesh(double Lx, double Ly, double Lz, int Nx, int Ny, int Nz, double sx, double sy,
                                double sz) {

    printf("Generating internal mesh...\n");


    // Auxiliary variables declaration
    double nodeX, nodeY, nodeZ;
    std::vector<int> vertices;
    Element auxElement;
    Face auxFace;
    Boundary auxBoundary;

    // Generation of the mesh vector of nodes
    for (int k = 0; k < Nz + 1; ++k) {
        for (int i = 0; i < Ny + 1; ++i) {
            for (int j = 0; j < Nx + 1; ++j) {

                // x component of the node
                if (sx == 0) {
                    nodeX = Lx/Nx*j;
                } else {
                    nodeX = 0 + 0.5*(1 + tanh(sx*(double(j)/Nx - 0.5))/tanh(sx/2))*(Lx - 0);
                }

                // y component of the node
                if (sy == 0) {
                    nodeY = Ly/Ny*i;
                } else {
                    nodeY = 0 + 0.5*(1 + tanh(sy*(double(i)/Ny - 0.5))/tanh(sy/2))*(Ly - 0);
                }

                // z component of the node
                if (sz == 0) {
                    nodeZ = Lz/Nz*k;
                } else {
                    nodeZ = 0 + 0.5*(1 + tanh(sz*(double(k)/Nz - 0.5))/tanh(sz/2))*(Lz - 0);
                }

                // Add the node to the mesh vector of nodes
                nodes.emplace_back(nodeX, nodeY, nodeZ);
                nNodes++;
            }
        }
    }


    // Generation of the elements of the mesh
    for (int k = 0; k < Nz; ++k) {
        for (int i = 0; i < Ny; ++i) {
            for (int j = 0; j < Nx; ++j) {

                // Computation of the vertices of each element
                vertices = iterators2ElementVertices(i, j, k, Nx, Ny, Nz);
                auxElement.iNodes = vertices;

                // Computation of the centroid of each element
                auxElement.centroid = computeHexahedronCentroid(nodes[vertices[0]], nodes[vertices[1]],
                                                                nodes[vertices[2]], nodes[vertices[3]],
                                                                nodes[vertices[4]], nodes[vertices[5]],
                                                                nodes[vertices[6]], nodes[vertices[7]]);

                // Computation of the volume of each element
                auxElement.Vf = computeHexahedronVolume(nodes[vertices[0]].x, nodes[vertices[1]].x,
                                                        nodes[vertices[0]].y, nodes[vertices[4]].y,
                                                        nodes[vertices[0]].z, nodes[vertices[3]].z);

                // Set the element type to a quad element
                auxElement.elementType = 12;

                // Add the interior element to the mesh vector of elements
                elements.push_back(auxElement);
                nInteriorElements++;
                nElements++;

            }
        }
    }


    // Generation of the yz faces (Perpendicular to the x-axis)
    for (int k = 0; k < Nz; ++k) {
        for (int i = 0; i < Ny; ++i) {
            for (int j = 0; j < Nx - 1; ++j) {

                // Computation of the vertices
                vertices = iterators2ElementVertices(i, j, k, Nx, Ny, Nz);
                auxFace.iNodes = {vertices[1], vertices[5], vertices[6], vertices[2]};

                // Computation of the owner and the neighbour element (face to element connectivity).
                // It is assumed that the owner is the element with the lowest index
                auxFace.iOwner = j + Nx*i + Nx*Ny*k;
                auxFace.iNeighbour = (j + 1) + Nx*i + Nx*Ny*k;

                // Element to face connectivity
                elements[auxFace.iOwner].iFaces.push_back(nFaces);
                elements[auxFace.iNeighbour].iFaces.push_back(nFaces);

                // Element to element connectivity
                elements[auxFace.iOwner].iElements.push_back(auxFace.iNeighbour);
                elements[auxFace.iNeighbour].iElements.push_back(auxFace.iOwner);

                // Computation of the surface and the surface vector. The surface vector points to the neighbour element
                auxFace.SfMag = computeRectangleArea(nodes[vertices[5]].y, nodes[vertices[1]].y,
                                                     nodes[vertices[2]].z, nodes[vertices[1]].z);
                auxFace.Sf = {auxFace.SfMag, 0, 0};

                // Add the auxiliary face to the mesh faces array
                faces.push_back(auxFace);
                nInteriorFaces++;
                nFaces++;
            }
        }
    }


    // Generation of the xz faces (Perpendicular to the y-axis)
    for (int k = 0; k < Nz; ++k) {
        for (int i = 0; i < Ny - 1; ++i) {
            for (int j = 0; j < Nx; ++j) {

                // Computation of the vertices
                vertices = iterators2ElementVertices(i, j, k, Nx, Ny, Nz);
                auxFace.iNodes = {vertices[4], vertices[7], vertices[6], vertices[5]};

                // Computation of the owner and the neighbour element (face to element connectivity).
                // It is assumed that the owner is the element with the lowest index
                auxFace.iOwner = j + Nx*i + Nx*Ny*k;
                auxFace.iNeighbour = j + Nx*(i + 1) + Nx*Ny*k;

                // Element to face connectivity
                elements[auxFace.iOwner].iFaces.push_back(nFaces);
                elements[auxFace.iNeighbour].iFaces.push_back(nFaces);

                // Element to element connectivity
                elements[auxFace.iOwner].iElements.push_back(auxFace.iNeighbour);
                elements[auxFace.iNeighbour].iElements.push_back(auxFace.iOwner);

                // Computation of the surface and the surface vector. The surface vector points to the neighbour element
                auxFace.SfMag = computeRectangleArea(nodes[vertices[7]].z, nodes[vertices[4]].z,
                                                     nodes[vertices[5]].x, nodes[vertices[4]].x);
                auxFace.Sf = {0, auxFace.SfMag, 0};

                // Add the auxiliary face to the mesh faces array
                faces.push_back(auxFace);
                nInteriorFaces++;
                nFaces++;
            }
        }
    }


    // Generation of the xy faces (Perpendicular to the z-axis)
    for (int k = 0; k < Nz - 1; ++k) {
        for (int i = 0; i < Ny; ++i) {
            for (int j = 0; j < Nx; ++j) {

                // Computation of the vertices
                vertices = iterators2ElementVertices(i, j, k, Nx, Ny, Nz);
                auxFace.iNodes = {vertices[3], vertices[2], vertices[6], vertices[7]};

                // Computation of the owner and the neighbour element (face to element connectivity).
                // It is assumed that the owner is the element with the lowest index
                auxFace.iOwner = j + Nx*i + Nx*Ny*k;
                auxFace.iNeighbour = j + Nx*i + Nx*Ny*(k + 1);

                // Element to face connectivity
                elements[auxFace.iOwner].iFaces.push_back(nFaces);
                elements[auxFace.iNeighbour].iFaces.push_back(nFaces);

                // Element to element connectivity
                elements[auxFace.iOwner].iElements.push_back(auxFace.iNeighbour);
                elements[auxFace.iNeighbour].iElements.push_back(auxFace.iOwner);

                // Computation of the surface and the surface vector. The surface vector points to the neighbour element
                auxFace.SfMag = computeRectangleArea(nodes[vertices[2]].x, nodes[vertices[3]].x,
                                                     nodes[vertices[7]].y, nodes[vertices[3]].y);
                auxFace.Sf = {0, 0, auxFace.SfMag};

                // Add the auxiliary face to the mesh faces array
                faces.push_back(auxFace);
                nInteriorFaces++;
                nFaces++;
            }
        }
    }


    // Loop over all the interior faces
    for (int i = 0; i < nInteriorFaces; ++i) {

        // Computation of the centroid
        faces[i].centroid = computeRectangleCentroid(nodes[faces[i].iNodes[0]], nodes[faces[i].iNodes[1]],
                                                     nodes[faces[i].iNodes[2]], nodes[faces[i].iNodes[3]]);

        // Computation of distances. Owner to neighbour & owner to face
        faces[i].dON = elements[faces[i].iNeighbour].centroid - elements[faces[i].iOwner].centroid;
        faces[i].dOf = faces[i].centroid - elements[faces[i].iOwner].centroid;
        faces[i].dONMag = faces[i].dON.mag();
        faces[i].dOfMag = faces[i].dOf.mag();

        // Computation of the linear interpolation factor to face
        faces[i].gf = (faces[i].dONMag - faces[i].dOfMag)/faces[i].dONMag;
    }


    printf("Internal mesh generated successfully!!\n");
}
