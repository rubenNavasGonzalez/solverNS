//
// Created by ruben on 19/11/23.
//

#ifndef FLOWBETWEENFLATPLATES_BOUNDARY_H
#define FLOWBETWEENFLATPLATES_BOUNDARY_H


class Boundary {
public:
    int startFace{};            // Starting face of the boundary
    int nBoundaryFaces{};       // Number of faces in the boundary


    // Boundary constructor and destructor
    Boundary();
    ~Boundary();
};


#endif //FLOWBETWEENFLATPLATES_BOUNDARY_H
