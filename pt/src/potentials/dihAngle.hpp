#pragma once
#include "abstract.hpp"

class DihAngle: public Potential {
public:
    DihAngle(const Config &_config,
             const State &_state,
             const PAtoms &_p_atoms);

    void init_and_register(VerList &verlet_list) override {}
    Vec3DArray calculate_forces(const VerList &verlet_list) override;
    bool is_enabled() const override;
private:
    bool enabled = false;
    const bool harmonic_version;
    const Scalar k_phi_harmonic, k_phi_1, k_phi_2;
    vector<Scalar> native_phi;
};
