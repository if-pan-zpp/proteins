#pragma once
#include "abstract.hpp"

class StructuredLJ: public Potential {
public:
    StructuredLJ(const Config &_config,
                 const State &_state,
                 const PAtoms &_p_atoms);

    Vec3DArray calculate_forces(const VerList &verlet_list) override;
    bool is_enabled() const override;
private:
    bool enabled = false;
    const Scalar force_cap = 1000.0;
    const Scalar r_cut, r_cut_sq; // Cut-off distance for this potential.

    struct native_contact {
        int i, j;
        Scalar sigma_sq;
        native_contact(int _i, int _j, Scalar _sigma_sq) :
            i(_i), j(_j), sigma_sq(_sigma_sq) {}
    };
    vector<native_contact> native_contacts;
};
