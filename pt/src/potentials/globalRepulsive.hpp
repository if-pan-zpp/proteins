#pragma once
#include "abstract.hpp"

class GlobalRepulsive: public Potential {
public:
    GlobalRepulsive(const Config &_config,
                    const State &_state,
                    const PAtoms &_p_atoms);

    void init_and_register(VerList &verlet_list) override;
    Vec3DArray calculate_forces(const VerList &verlet_list) override;
    bool is_enabled() const override;
private:
    bool enabled;
    Scalar cutoff;
    Scalar force_cap = 1000.0;
};
