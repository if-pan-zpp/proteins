#pragma once
#include "config.hpp"
#include "def/types.hpp"
using namespace std;

/*
 * Class that stores the state of the environment.
 * Examples: temperature, movement of walls, current step.
 *
 * We need to figure out which classes should be allowed to make changes here.
 */
class State {
public:
    State(const Config &config);
    bool is_running();
    void take_step();
    int max_steps = 0;
    int cur_step = 0;
    // only for prototype, generally temperature may change during simulation
    Scalar temperature = 0.35;
};
