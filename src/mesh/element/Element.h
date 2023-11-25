//
// Created by ruben on 19/11/23.
//

#ifndef FLOWBETWEENFLATPLATES_ELEMENT_H
#define FLOWBETWEENFLATPLATES_ELEMENT_H

#include <vector>
#include "../node/Node.h"


class Element {
public:
    std::vector<int> iNodes, iFaces;            // Tag of the nodes and faces forming the element
    std::vector<int> iElements;                 // Tag of the neighbour elements
    Node centroid;                              // Element centroid
    double Vf{};                                // Element volume
    int elementType{};                          // Type of element


    // Element constructor and destructor
    Element();
    ~Element();
};


#endif //FLOWBETWEENFLATPLATES_ELEMENT_H
