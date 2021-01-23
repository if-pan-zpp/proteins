#pragma once
#include "abstract.hpp"

class BondAngle: public Potential {
public:
    BondAngle(const Config &_config,
              const State &_state,
              const PAtoms &_p_atoms);

    def::Vec3DArray calculate_forces(const VerList &verlet_list) override;
    bool is_enabled() const override;
private:
    bool enabled = false;
};
