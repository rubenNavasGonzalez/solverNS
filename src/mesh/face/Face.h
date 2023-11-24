//
// Created by ruben on 19/11/23.
//

#ifndef FLOWBETWEENFLATPLATES_FACE_H
#define FLOWBETWEENFLATPLATES_FACE_H

#include <vector>
#include "../node/Node.h"
#include "../../math/geometricVector/GeometricVector.h"


class Face {
public:
    std::vector<int> iNodes;                // Tag of the nodes forming the face
    Node centroid;                          // Face centroid
    GeometricVector Sf;                     // Surface vector (oriented to owner element)
    double SfMag{};                         // Face area
    int iOwner{}, iNeighbour{};             // Tag of the owner and neighbour elements of the face
    int iOwnerFar{}, iNeighbourFar{};       // Tag of the far owner and far neighbour elements of the face
    GeometricVector dON, dOf;               // Owner to neighbour element and Owner to face distance vectors
    double dONMag{}, dOfMag{};              // Magnitude of the vectors above
    double gf{};                            // Face linear iterpolation factor Phi_f = gf*Phi_C + (1 - gf)*Phi_E
    int iPeriodicFace{};                    // Tag of the periodic face (just used in boundary faces)


    // Face constructor and destructor
    Face();
    ~Face();
};


#endif //FLOWBETWEENFLATPLATES_FACE_H
