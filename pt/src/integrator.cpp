#include "integrator.hpp"

Integrator::Integrator(const Config &_config,
                       const State &_state,
                       const PAtoms &_p_atoms) : state(_state) {}

void Integrator::take_step(PAtoms &p_atoms, const def::Vec3DArray &forces) {
    def::Scalar temp = state.temperature;
}
