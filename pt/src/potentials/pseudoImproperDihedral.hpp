#pragma once
#include "abstract.hpp"

class PseudoImproperDihedral: public Potential {
public:
    PseudoImproperDihedral(const Config &_config,
                    const State &_state,
                    const PAtoms &_p_atoms);

    void init_and_register(VerList &verlet_list) override;
    Vec3DArray calculate_forces(const VerList &verlet_list) override;
    bool is_enabled() const override;
private:
    bool enabled;
    Scalar cutoff;
    bool use_mj_matrix;
    const map<pair<string,string>, Scalar> &mj_matrix;
    const map<pair<string,string>, Scalar> &sigma_ss;
    bool pid_cos;
    Scalar alpha_bb_pos;
    Scalar alpha_bb_neg;
    Scalar alpha_ss;
    Scalar psi0_bb_pos;
    Scalar psi0_bb_neg;
    Scalar psi0_ss;
    Scalar rmin_pos;
    Scalar rmin_neg;
    Scalar contact_mltp;
    bool sink_pot;
    Scalar eps_bb;
    bool pid_barrier;

    Scalar min_lambda = 0.00005;
    Scalar min_norm = 0.01;
};
