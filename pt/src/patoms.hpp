#pragma once
#include "def/types.hpp"
#include "conf/constants.hpp"
#include "config.hpp"
#include <vector>
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
    vector<FatBool> connected; // 0 if (i, i+1) are in different chains, != 0 otherwise
    vector<u_int16_t> neighbours; // number of neighbours of residue
    Vec3DArray native_pos;
    vector<Scalar> native_distances;
    vector<pair<int, int>> native_contacts;

    // For testing, temporary. It's the number of times
    // ran2 was called before running the simulation.
    int rng_calls_cnt = 0; 
};
