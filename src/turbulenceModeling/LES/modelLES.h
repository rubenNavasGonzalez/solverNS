//
// Created by ruben on 9/02/24.
//

#ifndef SOLVERNS_MODELLES_H
#define SOLVERNS_MODELLES_H


// LES turbulence models enumeration
enum modelLES {
    None,
    Smagorinsky,
    WALE,
    Vreman,
    Verstappen,
    S3PQ,
    S3PR,
    S3QR
};


// Reference of the turbulence models
#include "modelSmagorinsky/modelSmagorinsky.h"
#include "modelWALE/modelWALE.h"
#include "modelVreman/modelVreman.h"
#include "modelVerstappen/modelVerstappen.h"
#include "modelS3PQR/modelS3PQR.h"

#endif //SOLVERNS_MODELLES_H