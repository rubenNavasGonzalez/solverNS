//
// Created by ruben on 15/12/23.
//

#include <cmath>
#include "computeTimeStepOrthogonal.h"
#include "../../fields/scalarField/ScalarField.h"


double computeTimeStepOrthogonal(const PolyMesh& theMesh, const VectorField& u, const ScalarField& nu, double f) {

    // Auxiliary variables
    int iOwner, iNeighbour;
    GeometricVector uElement;
    double SfMag, Vf, DeltaX, uMod;


    // Preallocate the time-step fields
    ScalarField DeltaT_c, DeltaT_d;


    // Loop over all the interior faces
    for (int i = 0; i < theMesh.nInteriorFaces; ++i) {

        iOwner = theMesh.faces[i].iOwner;
        iNeighbour = theMesh.faces[i].iNeighbour;

        SfMag = theMesh.faces[i].SfMag;
        Vf = theMesh.elements[iOwner].Vf;
        DeltaX = Vf/SfMag;
        uElement = u[iOwner];
        uMod = uElement.mag();

        DeltaT_c.push_back( 0.35*DeltaX/uMod );
        DeltaT_d.push_back( 0.2*pow(DeltaX,2)/nu[iOwner] );

        Vf = theMesh.elements[iNeighbour].Vf;
        DeltaX = Vf/SfMag;
        uElement = u[iNeighbour];
        uMod = uElement.mag();

        DeltaT_c.push_back( 0.35*DeltaX/uMod );
        DeltaT_d.push_back( 0.2*pow(DeltaX,2)/nu[iNeighbour] );
    }


    return f*std::min(DeltaT_c.min(), DeltaT_d.min());
}


double computeTimeStepOrthogonal(const PolyMesh& theMesh, const VectorField& u, double Pr, double f) {

    // Auxiliary variables
    int iOwner, iNeighbour;
    GeometricVector uElement;
    double SfMag, Vf, DeltaX, uMod;


    // Preallocate the time-step fields
    ScalarField DeltaT_c, DeltaT_d, DeltaT_t;


    // Loop over all the interior faces
    for (int i = 0; i < theMesh.nInteriorFaces; ++i) {

        iOwner = theMesh.faces[i].iOwner;
        iNeighbour = theMesh.faces[i].iNeighbour;

        SfMag = theMesh.faces[i].SfMag;
        Vf = theMesh.elements[iOwner].Vf;
        DeltaX = Vf/SfMag;
        uElement = u[iOwner];
        uMod = uElement.mag();

        DeltaT_c.push_back( 0.35*DeltaX/uMod );
        DeltaT_d.push_back( 0.2*pow(DeltaX,2)/Pr );
        DeltaT_t.push_back( 0.2*pow(DeltaX,2) );

        Vf = theMesh.elements[iNeighbour].Vf;
        DeltaX = Vf/SfMag;
        uElement = u[iNeighbour];
        uMod = uElement.mag();

        DeltaT_c.push_back( 0.35*DeltaX/uMod );
        DeltaT_d.push_back( 0.2*pow(DeltaX,2)/Pr );
        DeltaT_t.push_back( 0.2*pow(DeltaX,2) );
    }


    return f*std::min(DeltaT_c.min(), std::min(DeltaT_d.min(), DeltaT_t.min()));
}