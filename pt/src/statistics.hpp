#pragma once
#include "config.hpp"
#include "state.hpp"
#include "patoms.hpp"
using namespace std;

class Statistics {
public:
    Statistics(const Config &_config,
               const State &_state,
               const PAtoms &_p_atoms);

    void take_step();

private:
    const State &state;
    const PAtoms &p_atoms;
};
