#pragma once
#include "abstract.hpp"

class LocalRepulsive: public Potential {
public:
    LocalRepulsive(const Config &_config,
             const State &_state,
             const PAtoms &_p_atoms);

    Vec3DArray calculate_forces(const VerList &verlet_list) override;
    bool is_enabled() const override;
private:
    bool enabled;
    Scalar repulsive_cutoff;
};