//
// Created by ruben on 19/11/23.
//

#ifndef FLOWBETWEENFLATPLATES_POLYMESH_H
#define FLOWBETWEENFLATPLATES_POLYMESH_H

#include <string>
#include "../face/Face.h"
#include "../element/Element.h"
#include "../boundary/Boundary.h"


class PolyMesh {
public:
    int nNodes{};                                                   // Number of nodes
    int nInteriorElements{}, nBoundaryElements{}, nElements{};      // Number of elements
    int nInteriorFaces{}, nBoundaryFaces{}, nFaces{};               // Number of faces
    int nBoundaries{};                                              // Number of boundaries
    std::vector<Node> nodes;                                        // Nodes of the mesh
    std::vector<Face> faces;                                        // Faces of the mesh
    std::vector<Element> elements;                                  // Elements of the mesh
    std::vector<Boundary> boundaries;                               // Boundaries of the mesh


    // PolyMesh constructor and destructor
    PolyMesh();
    ~PolyMesh();


    // PolyMesh methods
    // Generate internal mesh
    void generatePolyMesh(double Lx, double Ly, double Lz, int Nx, int Ny, int Nz, double sx, double sy, double sz);
    // Generate boundary mesh (it must be used after the internal mesh generation)
    void generateBoundaryMesh(int Nx, int Ny, int Nz);
    // Write the mesh to a .VTK file
    void writeMesh2VTK(const std::string& filename) const;
};


#endif //FLOWBETWEENFLATPLATES_POLYMESH_H
