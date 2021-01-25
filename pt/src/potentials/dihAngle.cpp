#include "dihAngle.hpp"
#include <Eigen/Dense>
#include <iostream>

DihAngle::DihAngle(const Config &_config,
                   const State &_state,
                   const PAtoms &_p_atoms) :
    Potential(_config, _state, _p_atoms),
    k_phi_harmonic(_config.k_phi_harmonic),
    k_phi_1(_config.k_phi_1),
    k_phi_2(_config.k_phi_2),
    harmonic_version(_config.harmonic_dih_pot) {

    enabled = _config.dih_angle_pot && _p_atoms.n > 3;

    native_phi.resize(_p_atoms.n - 3, 0.0);
    const Vec3DArray &pos = p_atoms.native_pos;
    for (size_t i = 0; i + 3 < (size_t) p_atoms.n; ++i) {
        if (p_atoms.connected[i + 0] &&
            p_atoms.connected[i + 1] &&
            p_atoms.connected[i + 2]) {

            Vec3D v0 = pos.row(i + 1) - pos.row(i + 0);
            Vec3D v1 = pos.row(i + 2) - pos.row(i + 1);
            Vec3D v2 = pos.row(i + 3) - pos.row(i + 2);

            Vec3D normal_012 = v0.cross(v1).normalized();
            Vec3D normal_123 = v1.cross(v2).normalized();

            Scalar phi = acos(normal_012.dot(normal_123));
            if (normal_012.dot(v2) < 0.) phi = -phi;

            native_phi[i] = phi;
        }
    }
}

bool DihAngle::is_enabled() const {
    return enabled;
}

Vec3DArray DihAngle::calculate_forces(const VerList &verlet_list) {
    Vec3DArray forces = Vec3DArray::Zero(p_atoms.n, 3);
    const Vec3DArray &pos = p_atoms.der[0];

    Scalar energy = 0.0;
    Eigen::Array<Scalar, 4, 3> d_phi_d_pos; 

    for (size_t i = 0; i + 3 < (size_t) p_atoms.n; ++i) {
        if (p_atoms.connected[i + 0] &&
            p_atoms.connected[i + 1] &&
            p_atoms.connected[i + 2]) {

            Vec3D v0 = pos.row(i + 1) - pos.row(i + 0);
            Vec3D v1 = pos.row(i + 2) - pos.row(i + 1);
            Vec3D v2 = pos.row(i + 3) - pos.row(i + 2);

            Vec3D normal_012 = v0.cross(v1);
            Vec3D normal_123 = v1.cross(v2);
            Scalar normal_012_len_sq = normal_012.squaredNorm();
            Scalar normal_123_len_sq = normal_123.squaredNorm();

            Scalar phi = acos(normal_012.dot(normal_123) /
                              sqrt(normal_012_len_sq * normal_123_len_sq));
            if (normal_012.dot(v2) < 0.) phi = -phi;

            Scalar v1_len_sq = v1.squaredNorm();
            Scalar v1_len = sqrt(v1_len_sq);

            d_phi_d_pos.row(0) = -normal_012 * v1_len / normal_012_len_sq;
            d_phi_d_pos.row(3) = normal_123 * v1_len / normal_123_len_sq;

            // TODO: figure out why this.
            Scalar aux1 = -v0.dot(v1);
            Scalar aux2 = -v2.dot(v1);

            Eigen::Array<Scalar, 1, 3> df;
            df = d_phi_d_pos.row(0) * aux1 - d_phi_d_pos.row(3) * aux2;
            df /= v1_len_sq;

            d_phi_d_pos.row(1) = -d_phi_d_pos.row(0) + df;
            d_phi_d_pos.row(2) = -d_phi_d_pos.row(3) - df;

            phi -= native_phi[i];
            Scalar d_V_d_phi;
            if (harmonic_version) {
                energy += 0.5 * k_phi_harmonic * phi * phi;
                // This minus sign is almost definitely a bug, but it's also in cg.f 
                // TODO: prove that it's a bug.
                d_V_d_phi = -k_phi_harmonic * phi;
            }
            else {
                energy += k_phi_1 * (1.0 - cos(phi)) + k_phi_2 * (1.0 - cos(3. * phi));
                d_V_d_phi = k_phi_1 * sin(phi) + k_phi_2 * 3. * sin(3. * phi);
            }

            forces.block<4,3>(i,0) -= d_V_d_phi * d_phi_d_pos;
        }
    }

    return forces;
}
