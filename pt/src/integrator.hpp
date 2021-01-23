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

    void take_step(PAtoms &p_atoms, const def::Vec3DArray &forces);

private:
    const State &state;
};