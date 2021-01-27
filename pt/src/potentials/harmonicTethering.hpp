#pragma once
#include "abstract.hpp"

class HarmonicTethering: public Potential {
public:
    HarmonicTethering(const Config &_config,
             const State &_state,
             const PAtoms &_p_atoms);

    Vec3DArray calculate_forces(const VerList &verlet_list) override;
    bool is_enabled() const override;
private:
    bool enabled;
    Scalar H1;
    Scalar H2;
    Scalar force_cap = 1000.0;
};
