//
// Created by ruben on 21/01/24.
//

#ifndef SOLVERNS_WRITETURBULENTCHANNELFLOWDATA2CSV_H
#define SOLVERNS_WRITETURBULENTCHANNELFLOWDATA2CSV_H

#include "../../../fields/scalarField/ScalarField.h"
#include "../../../fields/vectorField/VectorField.h"
#include "../../../fields/tensorField/TensorField.h"


void writeTurbulentChannelFlowData2CSV(const PolyMesh& theMesh, const ScalarField& p, const VectorField& u, const VectorField& omega,
                                       const ScalarField& nut, const TensorField& gradU, double t);

#endif //SOLVERNS_WRITETURBULENTCHANNELFLOWDATA2CSV_H
