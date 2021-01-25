#pragma once
#include "config.hpp"
#include "state.hpp"
#include "patoms.hpp"
#include "def/types.hpp"
using namespace std;

/*
 * This deals with langevin dynamics (damping and thermal noise)
 * as well as the predictor-corrector algorithm.
 */
class Integrator {
public:
    Integrator(const Config &_config,
               const State &_state,
               const PAtoms &_p_atoms);

    void take_step(PAtoms &p_atoms, Vec3DArray &forces);

private:
    const State &state;
    const Scalar gamma;
    const Scalar delta;
    const int first_index;
    const int residues;
    const array<Scalar, DER_ORDER + 1> pred_corr_params =
    {3./16,  251./360,  1.,  11./18,  1./6,  1./60};
    const array<array<Scalar, DER_ORDER + 1>, DER_ORDER + 1> newton_symbol = 
    {{{{1, 1, 1, 1, 1,  1}},
      {{0, 1, 2, 3, 4,  5}},
      {{0, 0, 1, 3, 6, 10}},
      {{0, 0, 0, 1, 4, 10}},
      {{0, 0, 0, 0, 1,  5}},
      {{0, 0, 0, 0, 0,  1}}}};
};
