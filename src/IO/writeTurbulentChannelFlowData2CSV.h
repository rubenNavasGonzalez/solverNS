//
// Created by ruben on 21/01/24.
//

#ifndef SOLVERNS_WRITETURBULENTCHANNELFLOWDATA2CSV_H
#define SOLVERNS_WRITETURBULENTCHANNELFLOWDATA2CSV_H

#include "../fields/scalarField/ScalarField.h"
#include "../fields/vectorField/VectorField.h"


void writeTurbulentChannelFlowData2CSV(const ScalarField& p, const VectorField& u, const VectorField& omega, double uBulk, double t);

#endif //SOLVERNS_WRITETURBULENTCHANNELFLOWDATA2CSV_H
