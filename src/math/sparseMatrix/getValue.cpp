//
// Created by ruben on 6/12/23.
//

#include "SparseMatrix.h"


std::vector<double> SparseMatrix::getValue(std::vector<int> ij) {

    if(ij[0] < ij[1]) {

        if( !upperIndex.empty() ) {

            for(int i = 0; i <= upperIndex.size()-1; i++) {

                if( upperIndex[i][0] == ij[0] && upperIndex[i][1] == ij[1] ) {

                    return {upperValue[i], double(i)};
                }
            }
        } else {

            return {0, -1};
        }
    } else if(ij[0] > ij[1]) {

        if( !lowerIndex.empty() ) {

            for(int i = 0; i <= lowerIndex.size()-1; i++) {

                if( lowerIndex[i][0] == ij[0] && lowerIndex[i][1] == ij[1] ) {

                    return {lowerValue[i], double(i)};
                }
            }
        } else {

            return {0, -1};
        }
    } else {

        if( !diagIndex.empty() ) {

            for(int i = 0; i <= diagIndex.size()-1; i++) {

                if( diagIndex[i][0] == ij[0] && diagIndex[i][1] == ij[1] ) {

                    return {diagValue[i], double(i)};
                }
            }
        } else {

            return {0, -1};
        }
    }

    return {0, -1};
}