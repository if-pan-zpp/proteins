#include "pseudoImproperDihedral.hpp"


PseudoImproperDihedral::PseudoImproperDihedral(const Config &_config,
                                            const State &_state,
                                            const PAtoms &_p_atoms) :
    Potential(_config, _state, _p_atoms), mj_matrix(_config.mj_matrix),
    sigma_ss(_config.sigma_ss), residue_types(_config.residue_types) {
        enabled = _config.pseudo_improper_dihedral_pot;
        cutoff = _config.all_potentials_r_cut;
        use_mj_matrix = _config.use_mj_matrix;
        pid_cos = _config.pid_cos;
        alpha_bb_pos = _config.alpha_bb_pos;
        alpha_bb_neg = _config.alpha_bb_pos;
        alpha_ss = _config.alpha_ss;
        psi0_bb_pos = _config.psi0_bb_pos;
        psi0_bb_neg = _config.psi0_bb_neg;
        psi0_ss = _config.psi0_ss;
        rmin_pos = _config.pid_rmin_pos;
        rmin_neg = _config.pid_rmin_neg;
        contact_mltp = _config.contact_mltp;
        sink_pot = _config.sink_pot;
        eps_bb = _config.eps_bb;
        pid_barrier = _config.pid_barrier;
        pid_electrostatics = _config.pid_electrostatics;
        elektr_screen = _config.elektr_screen;
        coul = _config.coul;
        ele_perm_const = _config.ele_perm_const;
    }


void PseudoImproperDihedral::init_and_register(VerList &verlet_list) {
    verlet_list.register_eps(cutoff);
}

Vec3DArray PseudoImproperDihedral::calculate_forces(
                                            const VerList &verlet_list){
    Vec3DArray forces = Vec3DArray::Zero(p_atoms.n, 3);
    const Vec3DArray &positions = p_atoms.der[0];
    const vector<string> &names = p_atoms.residue_name;

    pair<VerIt, VerIt> list = verlet_list.get_verlet_list(cutoff);

    for (VerIt it = list.first; it != list.second; it++) {
        int i = it -> first;
        int j = it -> second;

        bool norm_too_small = false;
        Scalar eps_mj;
        if(use_mj_matrix)
            eps_mj = mj_matrix.at({names[i], names[j]});
        else
            eps_mj = 1.;

        // TODO implement PBC distance
        Vec3D diff_vec = positions.row(j) - positions.row(i);
        Scalar dist = diff_vec.norm();
        Scalar sq_dist = dist * dist;

        // TODO update residues' neighbour count

        if(sq_dist > cutoff * cutoff)   continue;

        if(names[i] == "PRO" || names[j] == "PRO")  continue;

        Scalar ss_lambda[2];
        Scalar bb_lambda[2];
        Scalar alpha1[2];
        Scalar alpha2[2];
        Scalar alpha3[2]; // a(_,4) in cg.f
        Scalar r_min;
        Vec3D f_var[2][4]; // TODO find out what it represents and name it properly
        Scalar dvdp[2]; // TODO find descriptive name for it

        for(size_t nr = 0; nr <=1; nr++) {
            int i1 = i;
            int i2 = j;
            if(nr == 1) swap(i1,i2);
            // TODO check if i1-1, i1+1 etc. may be out of bounds
            Vec3D v1 = positions.row(i1) - positions.row(i1+1);
            Vec3D v2 = positions.row(i1-1) - positions.row(i1+1);
            Vec3D v3 = positions.row(i1-1) - positions.row(i2);
            Vec3D v4 = v1.cross(v2);
            Vec3D v5 = v2.cross(v3);
            Scalar v4_norm_sq = v4.dot(v4);
            Scalar v5_norm_sq = v5.dot(v5);

            if(v4_norm_sq < min_norm || v5_norm_sq < min_norm) {
                norm_too_small = true;
                continue;
            }

            Scalar cospsi = v5.dot(v4) / sqrt(v5_norm_sq * v4_norm_sq);
            Scalar psi = acos(cospsi);
            if(v1.dot(v5) < 0)  psi *= -1.;

            alpha1[nr] = alpha_bb_pos * (psi - psi0_bb_pos);
            if(alpha1[nr] < M_PI && alpha1[nr] > -M_PI) {
                if(pid_cos) {
                    ss_lambda[nr] = 0.5 * (cos(alpha1[nr]) + 1);
                } else {
                    Scalar alpha_sq = (alpha1[nr] / M_PI) * (alpha1[nr] / M_PI);
                    ss_lambda[nr] = 1. - alpha_sq 
                                    / (2. * alpha_sq -  2. * abs((alpha1[nr] / M_PI)) + 1.);
                }
            }

            if(abs(i1 - i2) == 3) {
                if(nr == 0) {
                    alpha2[nr] = alpha_bb_pos * (psi - psi0_bb_pos);
                    alpha3[nr] = alpha_bb_pos;
                    r_min = rmin_pos;
                } else {
                    alpha2[nr] = alpha_bb_neg * (psi - psi0_bb_neg);
                    alpha3[nr] = alpha_bb_neg;
                    r_min = rmin_neg;
                }
            } else {
                alpha2[nr] = alpha_bb_pos * (psi - psi0_bb_pos);
                alpha3[nr] = alpha_bb_pos;
                r_min = rmin_pos;
                if(alpha2[nr] > M_PI || alpha2[nr] < -M_PI) {
                    alpha2[nr] = alpha_bb_neg * (psi - psi0_bb_neg);
                    alpha3[nr] = alpha_bb_neg;
                    r_min = rmin_neg;
                }
            }

            if(alpha2[nr] < M_PI && alpha2[nr] > -M_PI) {
                if(pid_cos) {
                    bb_lambda[nr] = 0.5 * (cos(alpha2[nr]) + 1);
                } else {
                    Scalar alpha_sq = (alpha2[nr] / M_PI) * (alpha2[nr] / M_PI);
                    bb_lambda[nr] = 1. - alpha_sq 
                                    / (2. * alpha_sq -  2. * abs((alpha2[nr] / M_PI)) + 1.);
                }
            }

            if(bb_lambda[nr] > min_lambda || ss_lambda[nr] > min_lambda) {
                f_var[nr][0] = v4 * v2.norm() /  v4_norm_sq;
                f_var[nr][1] = - v5 * v2.norm() /  v5_norm_sq;
                Vec3D df = (f_var[nr][0] * v1.dot(v2) - f_var[nr][1] * v3.dot(v2)) / v2.squaredNorm();
                f_var[nr][2] = - f_var[nr][0] + df;
                f_var[nr][3] = - f_var[nr][1] - df;
            }
        }

        if(norm_too_small)  continue;

        Scalar ss_lambda_ = ss_lambda[0] * ss_lambda[1];
        Scalar bb_lambda_ = bb_lambda[0] * bb_lambda[1];

        if(ss_lambda_ < min_lambda && bb_lambda_ < min_lambda)    continue;

        Scalar force;
        Scalar lj_energy;

        if(bb_lambda_ > min_lambda) {
            if(dist < r_min * contact_mltp) {
                //TODO update global contact count
            }
            if(sink_pot && dist < r_min * pow(2., 1. / 6.)) {
                lj_energy = - eps_bb;
            }
            else if(pid_barrier) {
                Scalar rsi = r_min *  pow(2., 1. / 6.) / dist;
                Scalar r6 = pow(rsi, 6.);
                lj_energy = r6 * (4. * r6 - 18. * rsi + 13.) * eps_bb;
                force += 6. * r6 * (21. * rsi - 8. * r6 - 13.)
                                        / dist * bb_lambda_ * eps_bb;
            } else {
                Scalar rsi = r_min  / dist;
                Scalar r6 = pow(rsi, 6.);
                lj_energy = 4. * r6 * (1. - r6) * eps_bb;
                force += 24. * r6 * (1. - 2. * r6) / dist * bb_lambda_ * eps_bb;
            }
            //TODO update global potential energy

            for(size_t nr = 0; nr <=1; nr++) {
                size_t other_nr = 1 - nr;
                if(pid_cos) {
                    dvdp[nr] -= 0.5 * alpha3[nr] * sin(alpha2[nr]) * bb_lambda[other_nr] * lj_energy;
                } else {
                    if(alpha2[nr] > 0.) {
                        Scalar dgdx = 2. * alpha2[nr] * (alpha2[nr] - 1.);
                        Scalar dgdx2 = (dgdx + 1.) * (dgdx + 1.);
                        dvdp[nr] += alpha3[nr] * dgdx / dgdx2 * bb_lambda[other_nr] * lj_energy;
                    } else {
                        Scalar dgdx = 2. * alpha2[nr] * (alpha2[nr] + 1.);
                        Scalar dgdx2 = (dgdx - 1.) * (dgdx - 1.);
                        dvdp[nr] += alpha3[nr] * dgdx / dgdx2 * bb_lambda[other_nr] * lj_energy;
                    }
                }
            }
        }
        //TODO consider moving analogous pieces of code for bb and ss to external method

        if(eps_mj > min_lambda && ss_lambda_ > min_lambda) {
            r_min = abs(sigma_ss.at({names[i], names[j]}));
            if(dist < r_min * contact_mltp) {
                //TODO update global contact count
            }
            int i_type = residue_types.at(names[i]);
            int j_type = residue_types.at(names[j]);
            bool i_charged = i_type == polar_pos || i_type == polar_neg;
            bool j_charged = j_type == polar_pos || j_type == polar_neg;
            if(pid_electrostatics && i_charged && j_charged){
                Scalar expc = exp(-1. / elektr_screen);
                Scalar coulpotcoeff = expc * coul;
                if(i_type == j_type)    coulpotcoeff *= -1.;
                if(ele_perm_const) {
                    lj_energy = coulpotcoeff / dist;
                    force += lj_energy * (-1. / elektr_screen - 1.) / dist;
                } else {
                    lj_energy = coulpotcoeff / sq_dist;
                    force += lj_energy * (-1. / elektr_screen - 2.) / dist;
                }
            }
            else if(sink_pot && dist < r_min * pow(2., 1. / 6.)) {
                lj_energy = - eps_mj;
            }
            else if(pid_barrier) {
                Scalar rsi = r_min *  pow(2., 1. / 6.) / dist;
                Scalar r6 = pow(rsi, 6.);
                lj_energy = r6 * (4. * r6 - 18. * rsi + 13.) * eps_mj;
                force += 6. * r6 * (21. * rsi - 8. * r6 - 13.) 
                                        / dist * ss_lambda_ * eps_mj;
            } else {
                Scalar rsi = r_min  / dist;
                Scalar r6 = pow(rsi, 6.);
                lj_energy = 4. * r6 * (1. - r6) * eps_mj;
                force += 24. * r6 * (1. - 2. * r6) / dist * ss_lambda_ * eps_mj;
            }
            //TODO update global potential energy

            for(size_t nr = 0; nr <=1; nr++) {
                size_t other_nr = 1 - nr;
                if(pid_cos) {
                    dvdp[nr] -= 0.5 * alpha_ss* sin(alpha1[nr]) * ss_lambda[other_nr] * lj_energy;
                } else {
                    if(alpha1[nr] > 0.) {
                        Scalar dgdx = 2. * alpha1[nr] * (alpha1[nr] - 1.);
                        Scalar dgdx2 = (dgdx + 1.) * (dgdx + 1.);
                        dvdp[nr] += alpha_ss * dgdx / dgdx2 * ss_lambda[other_nr] * lj_energy;
                    } else {
                        Scalar dgdx = 2. * alpha1[nr] * (alpha1[nr] + 1.);
                        Scalar dgdx2 = (dgdx - 1.) * (dgdx - 1.);
                        dvdp[nr] += alpha_ss * dgdx / dgdx2 * ss_lambda[other_nr] * lj_energy;
                    }
                }
            }
        }

        force /= -dist;
        Vec3D rep = force * diff_vec;

        forces.row(i) += rep.array();
        forces.row(j) -= rep.array();

        for(size_t nr = 0; nr <=1; nr++) {
            int i1 = i;
            int i2 = j;
            if(nr == 1) swap(i1,i2);
            forces.row(i1) -= dvdp[nr] * f_var[nr][0].array();
            forces.row(i1+1) -= dvdp[nr] * f_var[nr][1].array();
            forces.row(i1-1) -= dvdp[nr] * f_var[nr][2].array();
            forces.row(i2) -= dvdp[nr] * f_var[nr][3].array();
        }
    }
    return forces;
}


bool PseudoImproperDihedral::is_enabled() const {
    return enabled;
}