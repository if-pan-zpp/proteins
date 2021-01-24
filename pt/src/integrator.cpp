#include "integrator.hpp"

Integrator::Integrator(const Config &_config,
                       const State &_state,
                       const PAtoms &_p_atoms) : state(_state) {}

void Integrator::take_step(PAtoms &p_atoms, const Vec3DArray &forces) {
    Scalar temp = state.temperature;
}
