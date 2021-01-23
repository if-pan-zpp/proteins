#include "state.hpp"

State::State(const Config &config) {
    max_steps = config.max_steps;
}
bool State::is_running() {
    return cur_step < max_steps;
}
void State::take_step() {
    cur_step++;
}
