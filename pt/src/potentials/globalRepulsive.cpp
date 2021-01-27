#include "globalRepulsive.hpp"
#include <Eigen/Dense>
#include <iostream>

GlobalRepulsive::GlobalRepulsive(const Config &_config,
                                 const State &_state,
                                 const PAtoms &_p_atoms) :
    Potential(_config, _state, _p_atoms) {

    enabled = true;
    cutoff = _config.repulsive_cutoff;
}

void GlobalRepulsive::init_and_register(VerList &verlet_list) {
    verlet_list.register_eps(cutoff);
}

Vec3DArray GlobalRepulsive::calculate_forces(const VerList &verlet_list) {
    Vec3DArray forces = Vec3DArray::Zero(p_atoms.n, 3);
    const Vec3DArray &positions = p_atoms.der[0];

    pair<VerIt, VerIt> list = verlet_list.get_verlet_list(cutoff);
    Scalar sigma = cutoff / pow(2., 1./6.);
    for (VerIt it = list.first; it != list.second; it++) {
        int i = it -> first;
        int j = it -> second;
        Vec3D diff_vec = positions.row(j) - positions.row(i);
        Scalar dist = diff_vec.norm();
        Scalar sq_dist = dist * dist;

        if(sq_dist < cutoff * cutoff) {
            Scalar rsi = sigma / dist;
            Scalar r6 = pow(rsi, 6.);
            Scalar energy = 4 * r6 * (r6 - 1) + 1;
            Scalar force = 24 * r6 * (1 - 2 * r6) / dist;
            if (force > force_cap) force = force_cap;
            if (force < -force_cap) force = -force_cap;
            force /= -dist;
            forces.row(i) -= diff_vec.array() * force;
            forces.row(j) += diff_vec.array() * force;
        }
    }

    return forces;
}

bool GlobalRepulsive::is_enabled() const {
    return enabled;
}
