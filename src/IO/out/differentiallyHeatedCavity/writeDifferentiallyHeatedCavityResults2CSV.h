//
// Created by ruben on 7/04/24.
//

#ifndef SOLVERNS_WRITEDIFFERENTIALLYHEATEDCAVITYRESULTS2CSV_H
#define SOLVERNS_WRITEDIFFERENTIALLYHEATEDCAVITYRESULTS2CSV_H

#include "../../../fields/vectorField/VectorField.h"
#include "../../../fields/scalarField/ScalarField.h"


void writeDifferentiallyHeatedCavityResults2CSV(const PolyMesh& theMesh, const VectorField& u, const ScalarField& p, const ScalarField& T,
                                                const ScalarField& nut, double t);

#endif //SOLVERNS_WRITEDIFFERENTIALLYHEATEDCAVITYRESULTS2CSV_H
