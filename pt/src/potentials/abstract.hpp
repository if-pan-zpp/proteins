#pragma once
#include "../config.hpp"
#include "../state.hpp"
#include "../patoms.hpp"
#include "../verlist.hpp"
#include "../statistics.hpp"
#include <memory>
using namespace std;

class Potential {
public:
    Potential(const Config &_config,
              const State &_state,
              const PAtoms &_p_atoms) :
        state(_state), p_atoms(_p_atoms) {}

    virtual void init_and_register(VerList &verlet_list) {}

    /*
     * We split one step into those three functions, because
     * in the future we may want to be able to write calculate_forces
     * in a way that can run concurrently with other calculate_forces calls
     * (from this and other potentials). But some elements need synchronization anyway.
     */
    virtual void init_step() {}
    virtual Vec3DArray calculate_forces(const VerList &verlet_list) = 0;
    virtual void finish_step(Statistics &statistics) {}

    virtual bool is_enabled() const = 0;

protected:
    const State &state;
    const PAtoms &p_atoms;
};

vector<unique_ptr<Potential>> all_supported_potentials(const Config &config,
                                                       const State &state,
                                                       const PAtoms &p_atoms);
