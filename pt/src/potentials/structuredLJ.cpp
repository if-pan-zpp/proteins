#include "structuredLJ.hpp"
#include <iostream>
#include <iomanip>

StructuredLJ::StructuredLJ(const Config &_config,
                           const State &_state,
                           const PAtoms &_p_atoms) :
    Potential(_config, _state, _p_atoms),
    r_cut(_config.all_potentials_r_cut),
    r_cut_sq(r_cut * r_cut) {

    enabled = _p_atoms.native_contacts.size() > 0;

    // Calculate sigmas squared for native contacts.
    const Vec3DArray &pos = _p_atoms.native_pos;
    const Scalar scaling_constant = pow(0.5, 1./3);
    for (const pair<int, int> &contact : _p_atoms.native_contacts) {
        int i = contact.first;
        int j = contact.second;
        Scalar dist_sq = (pos.row(i) - pos.row(j)).matrix().squaredNorm();
        Scalar sigma_sq = dist_sq * scaling_constant;
        native_contacts.emplace_back(i, j, sigma_sq);
    }
}

bool StructuredLJ::is_enabled() const {
    return enabled;
}

Vec3DArray StructuredLJ::calculate_forces(const VerList &verlet_list) {
    Vec3DArray forces = Vec3DArray::Zero(p_atoms.n, 3);
    const Vec3DArray &pos = p_atoms.der[0];

    Scalar energy = 0.;

    for (const native_contact &contact : native_contacts) {
        Vec3D v = pos.row(contact.j) - pos.row(contact.i);
        Scalar dist_sq = v.squaredNorm();
        if (dist_sq > r_cut_sq) continue;

        Scalar r = sqrt(dist_sq);

        Scalar sigma_by_r_6 = contact.sigma_sq / dist_sq; 
        sigma_by_r_6 = sigma_by_r_6 * sigma_by_r_6 * sigma_by_r_6;
        
        energy += 4.0 * sigma_by_r_6 * (sigma_by_r_6 - 1.0);
        Scalar force = 24.0 * sigma_by_r_6 * (1.0 - 2.0 * sigma_by_r_6) / r;

        // If not for this check, we could avoid calculating r (one sqrt call).
        if (force > force_cap) force = force_cap;
        if (force < -force_cap) force = -force_cap;

        force /= r;
        forces.row(contact.i) += v.array() * force;
        forces.row(contact.j) -= v.array() * force;
    }

    // Send energy to statistics
    
    return forces;
}
