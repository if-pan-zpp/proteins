#include "harmonicTethering.hpp"
#include <Eigen/Dense>

HarmonicTethering::HarmonicTethering(const Config &_config,
                     const State &_state,
                     const PAtoms &_p_atoms) :
    Potential(_config, _state, _p_atoms) {

    enabled = true;
    H1 = _config.H1;
    H2 = _config.H1;

}

Vec3DArray HarmonicTethering::calculate_forces(const VerList &verlet_list) {
    Vec3DArray forces = Vec3DArray::Zero(p_atoms.n, 3);
    const Vec3DArray &positions = p_atoms.der[0];
    const vector<FatBool> connected = p_atoms.connected;
    size_t residues = p_atoms.n;
    vector<Scalar> native_distances = p_atoms.native_distances;

    for(size_t i = 0; i < residues - 1; i++) {
        if(connected[i]) {
            Vec3D diff_vec = positions.row(i + 1) - positions.row(i);
            Scalar dist = diff_vec.norm();
            Scalar dist_change = dist - native_distances[i];
            Scalar sq_dist_change = dist_change * dist_change;
            Scalar energy = (H1 + H2 * sq_dist_change) * sq_dist_change;
            Scalar force = (2 * H1 + 4 * H2 * dist_change) * dist_change;
            if(force > force_cap) {
                force = force_cap;
            } else if(force < -force_cap) {
                force = -force_cap;
            }
            force /= -dist;
            forces.row(i) += diff_vec.array() * force;
            forces.row(i+1) -= diff_vec.array() * force;
        }

    }
}

bool HarmonicTethering::is_enabled() const {
    return enabled;
}
