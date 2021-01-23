#include "dihAngle.hpp"

DihAngle::DihAngle(const Config &_config,
                   const State &_state,
                   const PAtoms &_p_atoms) :
    Potential(_config, _state, _p_atoms) {

    enabled = _config.dih_angle_pot;
}

void DihAngle::init_and_register(VerList &verlet_list) {
    // Right now, if this potential is enabled, this will fail,
    // because maximum eps is 10.0.
    verlet_list.register_eps(15.0);
}
def::Vec3DArray DihAngle::calculate_forces(const VerList &verlet_list) {
    return def::Vec3DArray::Zero(p_atoms.n, 3);
}
bool DihAngle::is_enabled() const {
    return enabled;
}
