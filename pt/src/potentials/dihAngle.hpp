#pragma once
#include "abstract.hpp"

class DihAngle: public Potential {
public:
    DihAngle(const Config &_config,
             const State &_state,
             const PAtoms &_p_atoms);

    void init_and_register(VerList &verlet_list) override;
    def::Vec3DArray calculate_forces(const VerList &verlet_list) override;
    bool is_enabled() const override;
private:
    bool enabled = false;
};
