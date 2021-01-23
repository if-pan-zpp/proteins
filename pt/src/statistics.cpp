#include "statistics.hpp"
#include <iostream>

Statistics::Statistics(const Config &_config,
                       const State &_state,
                       const PAtoms &_p_atoms)
    : state(_state), p_atoms(_p_atoms) {};

void Statistics::take_step() {
    cout << "We're in step " << state.cur_step << endl;

}
