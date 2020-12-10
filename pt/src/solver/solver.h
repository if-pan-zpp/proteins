#pragma once
#include "state/state.h"

namespace solver {
    class Solver {
    public:
        virtual void advance(state::State& state) = 0;
    };
}
