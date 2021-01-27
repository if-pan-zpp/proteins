#include "localRepulsive.hpp"
#include <Eigen/Dense>

LocalRepulsive::LocalRepulsive(const Config &_config,
                     const State &_state,
                     const PAtoms &_p_atoms) :
    Potential(_config, _state, _p_atoms) {

    enabled = true;
    repulsive_cutoff = _config.repulsive_cutoff;

}

Vec3DArray LocalRepulsive::calculate_forces(const VerList &verlet_list) {
    Vec3DArray forces = Vec3DArray::Zero(p_atoms.n, 3);
    const Vec3DArray &positions = p_atoms.der[0];
    const vector<FatBool> connected = p_atoms.connected;
    size_t residues = p_atoms.n;
    vector<Scalar> native_distances = p_atoms.native_distances;

    for(size_t i = 0; i < residues - 2; i++) {
        if(connected[i] & connected[i+1]) {
            Vec3D diff_vec = positions.row(i + 2) - positions.row(i);
            Scalar dist = diff_vec.norm();
            Scalar sq_dist = dist * dist;

            if(sq_dist < repulsive_cutoff * repulsive_cutoff) {
                Scalar rsi = repulsive_cutoff / ((pow(2., 1./6.)) * dist);
                Scalar r6 = pow(rsi, 6.);
                Scalar energy = 4 * r6 * (r6 - 1) + 1;
                Scalar force = 24 * r6 * (1 - 2 * r6) / dist;
                if (force > force_cap) force = force_cap;
                if (force < -force_cap) force = -force_cap;
                force /= -dist;
                forces.row(i) -= diff_vec.array() * force;
                forces.row(i+2) += diff_vec.array() * force;
            }
        }

    }

    return forces;
}

bool LocalRepulsive::is_enabled() const {
    return enabled;
}
