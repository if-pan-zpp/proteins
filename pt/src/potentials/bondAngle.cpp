#include "bondAngle.hpp"

BondAngle::BondAngle(const Config &_config,
                     const State &_state,
                     const PAtoms &_p_atoms) :
    Potential(_config, _state, _p_atoms) {

    enabled = _config.angle_pot;
}

Vec3DArray BondAngle::calculate_forces(const VerList &verlet_list) {
    return Vec3DArray::Zero(p_atoms.n, 3);
}
bool BondAngle::is_enabled() const {
    return enabled;
}
