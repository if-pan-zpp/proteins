#pragma once
#include "state/def.hpp"
#include "conf/constants.hpp"
#include "config.hpp"
using namespace std;

/*
 * Class for storing pseudo-atoms. It holds:
 * - positions and their derivatives
 * - masses, amino acid sequence etc.
 * - division into chains
 * - division into structured and unstructured parts
 * - method for calculating the distance between two pseudatoms
 *   (because it depends on periodic boundary conditions and is used a lot)
 * - ...
 *
 * After construction, only Integrator can make changes here.
 */
class PAtoms {
public:
    PAtoms(const Config &config);
    int n;
    array<Vec3DArray, DER_ORDER + 1> der; // positions' derivatives.

};
