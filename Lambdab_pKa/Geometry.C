#pragma once
#include <cmath>
#include <random>

#include "Rtypes.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"

// calculating the distance from the origin
Float_t Length(Float_t X, Float_t Y, Float_t Z) {
    return pow(X * X + Y * Y + Z * Z, 0.5);
}


// angle between two vectors
Float_t Angle(Float_t X1, Float_t Y1, Float_t Z1, Float_t X2, Float_t Y2, Float_t Z2) {
    Float_t angle = acos((X1 * X2 + Y1 * Y2 + Z1 * Z2) / (Length(X1, Y1, Z1) * Length(X2, Y2, Z2)));
    return angle;
}
