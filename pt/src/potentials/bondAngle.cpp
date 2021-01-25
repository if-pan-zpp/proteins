#include "bondAngle.hpp"
#include <Eigen/Dense>

#include <iostream>

BondAngle::BondAngle(const Config &_config,
                     const State &_state,
                     const PAtoms &_p_atoms) :
    Potential(_config, _state, _p_atoms) {

    enabled = _config.angle_pot && (p_atoms.n > 2);
    k_theta = _config.k_theta;

    const Vec3DArray &pos = p_atoms.native_pos;
    
    native_theta.resize(p_atoms.n - 2, 0.0);
    for (size_t i = 0; i < (size_t) (p_atoms.n - 2); ++i) {
        if (p_atoms.connected[i] && p_atoms.connected[i + 1]) {
            Vec3D v0 = pos.row(i + 1) - pos.row(i);
            Vec3D v1 = pos.row(i + 2) - pos.row(i + 1);

            Scalar v0_norm = v0.norm();
            Scalar v1_norm = v1.norm();
            Vec3D v0_x_v1 = v0.cross(v1);

            native_theta[i] = acos(-v0.dot(v1) / (v0_norm * v1_norm));
        }
    }
}

Vec3DArray BondAngle::calculate_forces(const VerList &verlet_list) {
    Vec3DArray forces = Vec3DArray::Zero(p_atoms.n, 3);
    const Vec3DArray &pos = p_atoms.der[0];

    Scalar energy = 0.0;
    Eigen::Array33d d_theta_d_pos;

    for (size_t i = 0; i < (size_t) (p_atoms.n - 2); ++i) {
        if (p_atoms.connected[i] && p_atoms.connected[i + 1]) {
            Vec3D v0 = pos.row(i + 1) - pos.row(i);
            Vec3D v1 = pos.row(i + 2) - pos.row(i + 1);

            Scalar v0_norm = v0.norm();
            Scalar v1_norm = v1.norm();
            Vec3D v0_x_v1 = v0.cross(v1);

            Scalar theta = acos(-v0.dot(v1) / (v0_norm * v1_norm));

            Vec3D grad0 = v0.cross(v0_x_v1).normalized() / v0_norm;
            Vec3D grad2 = v1.cross(v0_x_v1).normalized() / v1_norm;

            d_theta_d_pos.row(0) = grad0;
            d_theta_d_pos.row(1) = -grad0 - grad2;
            d_theta_d_pos.row(2) = grad2;

            // Ifs for different versions of this potential will be here.
            // This is the structured case:
            theta -= native_theta[i];
            energy += k_theta * theta * theta;

            Scalar d_V_d_theta = 2.0 * k_theta * theta;

            forces.block<3,3>(i,0) -= d_V_d_theta * d_theta_d_pos;
        }
    }

    // TODO: Send energy to statistics. 

    return forces;
}
bool BondAngle::is_enabled() const {
    return enabled;
}
