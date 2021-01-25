#include "integrator.hpp"
#include "util/ran2.hpp"
#include <iostream>

Integrator::Integrator(const Config &_config,
                       const State &_state,
                       const PAtoms &_p_atoms) :
    state(_state),
    gamma(_config.gamma),
    delta(_config.delta),
    first_index(_config.first_index),
    residues(_p_atoms.n) {

    for (int i = 0; i < _p_atoms.rng_calls_cnt; ++i) {
        util::ran2(); // Advance rng state so that the results match with cg.f.
    }
}

Scalar normalScalar() {
    Scalar r1 = util::ran2();
    Scalar r2 = util::ran2();
    return sqrt(-2.0 * log(r1)) * cos(2 * M_PI * r2);
}

void Integrator::take_step(PAtoms &p_atoms, Vec3DArray &forces) {
    Scalar temperature = state.temperature;

    // lang in cg.f
    const Scalar gamma2 = gamma / delta;
    const Scalar const2 = sqrt(2 * temperature * gamma * delta) * delta;

    Vec3DArray &velocities = p_atoms.der[1]; // pseudo-atoms' velocities
    const Vec3DArray &accelerations = p_atoms.der[2]; // pseudo-atoms' accelerations

    Vec3DArray noise = Vec3DArray::Zero(p_atoms.n, 3);
    for (size_t i = 0; i < 3 * p_atoms.n; ++i) {
        noise(i) = normalScalar();
    }
    velocities += const2 * noise;
    forces -= gamma2 * velocities;
    
    // corr in cg.f
    const Scalar deltsq = 0.5 * delta * delta;

    const Vec3DArray err = accelerations - deltsq * forces;
    for (int der_nr = 0; der_nr <= DER_ORDER; der_nr++) {
        p_atoms.der[der_nr] -= err * pred_corr_params[der_nr];
    }

    //predct in cg.f
    for(int der_nr = 0; der_nr <= DER_ORDER; der_nr++) {
        for(int next_der = der_nr + 1; next_der <= DER_ORDER; next_der++) {
            p_atoms.der[der_nr] +=
                newton_symbol[der_nr][next_der] * p_atoms.der[next_der];
        }
    }
}
